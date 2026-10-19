#include "DiagBase.H"
#include <AMReX.H>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>

namespace {
void
check(bool condition, const std::string& message)
{
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void
parameters(const std::string& name, amrex::Real position)
{
  amrex::ParmParse pp(name);
  pp.add("int", 2);
  pp.add("file", name + ".dat");
  pp.add("volume_fraction_field", "vfrac");
  pp.addarr("position", amrex::Vector<amrex::Real>(AMREX_SPACEDIM, position));
  pp.addarr("field_names", amrex::Vector<std::string>{"Temp", "pressure"});
}

struct Fixture
{
  amrex::Vector<amrex::Geometry> geoms;
  amrex::Vector<amrex::BoxArray> grids;
  amrex::Vector<amrex::DistributionMapping> dmaps;
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> state;
  amrex::Vector<std::string> vars{"pressure", "Temp", "vfrac"};

  Fixture()
  {
    amrex::Box domain(amrex::IntVect(4), amrex::IntVect(11));
    amrex::RealBox physical(
      {AMREX_D_DECL(1.0, 1.0, 1.0)}, {AMREX_D_DECL(3.0, 3.0, 3.0)});
    amrex::Array<int, AMREX_SPACEDIM> periodic{};
    for (int lev = 0; lev < 2; ++lev) {
      geoms.emplace_back(domain, &physical, 0, periodic.data());
      domain.refine(2);
    }
    grids.emplace_back(geoms[0].Domain());
    grids[0].maxSize(4);
    grids.emplace_back(amrex::Box(amrex::IntVect(16), amrex::IntVect(23)));
    grids[1].maxSize(4);
    rebuild();
  }

  void rebuild()
  {
    state.clear();
    dmaps.clear();
    for (int lev = 0; lev < 2; ++lev) {
      amrex::Vector<int> owners(grids[lev].size());
      for (int g = 0; g < grids[lev].size(); ++g) {
        owners[g] = (g + 1) % amrex::ParallelDescriptor::NProcs();
      }
      dmaps.emplace_back(owners);
      state.push_back(
        std::make_unique<amrex::MultiFab>(grids[lev], dmaps[lev], 3, 1));
      state.back()->setVal(-999.0);
      state.back()->setVal(lev == 0 ? 10.0 : 20.0, 0, 1, 0);
      state.back()->setVal(lev == 0 ? 100.0 : 200.0, 1, 1, 0);
      state.back()->setVal(1.0, 2, 1, 0);
    }
  }

  void sample(DiagBase& diag, int step, amrex::Real time)
  {
    diag.prepare(2, geoms, grids, dmaps, vars);
    diag.processDiag(step, time, amrex::GetVecOfConstPtrs(state), vars);
  }
};

amrex::Vector<amrex::Vector<amrex::Real>>
rows(const std::string& file)
{
  std::ifstream in(file);
  check(in.good(), "missing point history: " + file);
  amrex::Vector<amrex::Vector<amrex::Real>> result;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream record(line);
    amrex::Vector<amrex::Real> row;
    amrex::Real value;
    while (record >> value) {
      row.push_back(value);
    }
    check(record.eof(), "nonnumeric history record");
    result.push_back(row);
  }
  return result;
}

void
history()
{
  parameters("probe", 2.13);
  Fixture fixture;
  auto diag = DiagBase::create("DiagPoint");
  diag->init("probe", "probe");
  amrex::Vector<std::string> requested;
  diag->addVars(requested);
  check(
    std::find(requested.begin(), requested.end(), "pressure") !=
      requested.end(),
    "pressure not requested from PeleC");
  check(!diag->doDiag(0.1, 1) && diag->doDiag(0.2, 2), "step interval");
  fixture.sample(*diag, 2, 0.2);

  // Move the fine patch away from the point: the next record must use coarse
  // data.
  fixture.grids[1] =
    amrex::BoxArray(amrex::Box(amrex::IntVect(8), amrex::IntVect(15)));
  fixture.grids[1].maxSize(4);
  fixture.rebuild();
  fixture.sample(*diag, 4, 0.4);
  diag->close();
  diag.reset();

  // A new diagnostic instance simulates reopening the same history after
  // restart.
  diag = DiagBase::create("DiagPoint");
  diag->init("probe", "probe");
  fixture.sample(*diag, 4, 0.4);
  fixture.sample(*diag, 6, 0.6);
  diag->close();
  if (amrex::ParallelDescriptor::IOProcessor()) {
    auto data = rows("probe.dat");
    check(data.size() == 3, "restart must append without duplicate rows");
    const int level_col = 2 + AMREX_SPACEDIM;
    const int field_col = 3 + 2 * AMREX_SPACEDIM;
    for (int r = 0; r < 3; ++r) {
      check(data[r].size() == field_col + 2, "wrong output column count");
      check(data[r][0] == 2 * (r + 1), "wrong step");
      check(std::abs(data[r][1] - 0.2 * (r + 1)) < 1.e-12, "wrong time");
      check(
        data[r][level_col] == (r == 0 ? 1 : 0), "wrong finest covering level");
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        check(data[r][2 + d] == 2.13, "requested position changed");
        check(
          data[r][level_col + 1 + d] == (r == 0 ? 2.1875 : 2.125),
          "wrong containing-cell center");
      }
      check(data[r][field_col] == (r == 0 ? 200.0 : 100.0), "wrong Temp/order");
      check(
        data[r][field_col + 1] == (r == 0 ? 20.0 : 10.0),
        "wrong pressure/order");
    }
  }
}

void
boundary(const std::string& scenario)
{
  const amrex::Real position = scenario == "upper_boundary"
                                 ? 3.0
                                 : scenario == "lower_boundary" ? 1.0 : 2.0;
  parameters("probe", position);
  Fixture fixture;
  auto diag = DiagBase::create("DiagPoint");
  diag->init("probe", "probe");
  fixture.sample(*diag, 2, 0.2);
  if (amrex::ParallelDescriptor::IOProcessor()) {
    auto data = rows("probe.dat");
    check(data.size() == 1, "boundary record count");
    const int level_col = 2 + AMREX_SPACEDIM;
    const amrex::Real center =
      position == 3.0 ? 2.9375 : position == 1.0 ? 1.125 : 2.0625;
    check(data[0][level_col] == (position == 1.0 ? 0 : 1), "boundary level");
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      check(data[0][level_col + 1 + d] == center, "boundary cell center");
    }
  }
}

void
noVolumeFraction()
{
  parameters("probe", 2.13);
  amrex::ParmParse("probe").remove("volume_fraction_field");
  Fixture fixture;
  fixture.vars.pop_back();
  auto diag = DiagBase::create("DiagPoint");
  diag->init("probe", "probe");
  amrex::Vector<std::string> requested;
  diag->addVars(requested);
  check(requested.size() == 2, "unspecified volume fraction must not be required");
  fixture.sample(*diag, 2, 0.2);
  if (amrex::ParallelDescriptor::IOProcessor()) {
    auto data = rows("probe.dat");
    check(data.size() == 1 && data[0].back() == 20.0,
          "sampling without a volume fraction field");
  }
}

std::string
expectedError(const std::string& scenario)
{
  if (scenario == "outside")
    return "outside the physical domain";
  if (scenario == "bad_dimension")
    return "position must have";
  if (scenario == "nonfinite")
    return "position must be finite";
  if (scenario == "bad_interval")
    return "positive int";
  if (scenario == "period")
    return "per is not supported";
  if (scenario == "filters")
    return "filters are not supported";
  if (scenario == "empty_fields")
    return "field_names must not be empty";
  if (scenario == "missing_field")
    return "wasn't found";
  if (scenario == "covered")
    return "fully covered";
  if (scenario == "changed_fields")
    return "columns do not match";
  if (scenario == "changed_position")
    return "position does not match";
  if (scenario == "older_restart")
    return "older than the last record";
  if (scenario == "inconsistent_step")
    return "inconsistent step/time";
  if (scenario == "broken_history")
    return "invalid history record";
  if (scenario == "unwritable")
    return "cannot open";
  return "";
}

void
rejected(const std::string& scenario)
{
  parameters("probe", scenario == "outside" ? 3.01 : 2.13);
  amrex::ParmParse pp("probe");
  if (scenario == "bad_dimension") {
    pp.addarr("position", amrex::Vector<amrex::Real>{2.13});
  } else if (scenario == "nonfinite") {
    pp.addarr(
      "position",
      amrex::Vector<amrex::Real>(
        AMREX_SPACEDIM, std::numeric_limits<amrex::Real>::infinity()));
  } else if (scenario == "bad_interval") {
    pp.add("int", 0);
    pp.add("per", 0.1);
  } else if (scenario == "period") {
    pp.add("per", 0.1);
  } else if (scenario == "filters") {
    pp.addarr("filters", amrex::Vector<std::string>{"hot"});
    amrex::ParmParse filter("probe.hot");
    filter.add("field_name", "Temp");
    filter.add("value_greater", 0.0);
  } else if (scenario == "empty_fields") {
    pp.remove("field_names");
  } else if (scenario == "missing_field") {
    pp.addarr("field_names", amrex::Vector<std::string>{"missing"});
  } else if (scenario == "unwritable") {
    pp.add("file", "missing-directory/probe.dat");
  }
  Fixture fixture;
  auto diag = DiagBase::create("DiagPoint");
  diag->init("probe", "probe");
  if (scenario == "covered")
    fixture.state[1]->setVal(0.0, 2, 1, 0);
  fixture.sample(*diag, 2, 0.2);
  if (scenario == "older_restart" || scenario == "inconsistent_step") {
    diag.reset();
    diag = DiagBase::create("DiagPoint");
    diag->init("probe", "probe");
    fixture.sample(
      *diag, scenario == "older_restart" ? 1 : 3,
      scenario == "older_restart" ? 0.1 : 0.2);
  } else if (scenario == "changed_fields" || scenario == "changed_position") {
    diag.reset();
    if (scenario == "changed_fields") {
      pp.addarr("field_names", amrex::Vector<std::string>{"pressure", "Temp"});
    } else {
      pp.addarr("position", amrex::Vector<amrex::Real>(AMREX_SPACEDIM, 2.14));
    }
    diag = DiagBase::create("DiagPoint");
    diag->init("probe", "probe");
    fixture.sample(*diag, 4, 0.4);
  } else if (scenario == "broken_history") {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      std::ofstream out("probe.dat", std::ios::app);
      out << "4 0.4 2.13\n";
    }
    amrex::ParallelDescriptor::Barrier();
    diag.reset();
    diag = DiagBase::create("DiagPoint");
    diag->init("probe", "probe");
    fixture.sample(*diag, 6, 0.6);
  }
}
} // namespace

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  int status = 0;
  int caught = 0;
  std::string scenario;
  amrex::ParmParse().get("scenario", scenario);
  const auto expected = expectedError(scenario);
  if (amrex::ParallelDescriptor::IOProcessor())
    std::remove("probe.dat");
  amrex::ParallelDescriptor::Barrier();
  try {
    if (scenario == "no_volume_fraction")
      noVolumeFraction();
    else if (scenario == "history")
      history();
    else if (scenario.find("boundary") != std::string::npos)
      boundary(scenario);
    else
      rejected(scenario);
  } catch (const std::exception& error) {
    caught = 1;
    if (
      expected.empty() ||
      std::string(error.what()).find(expected) == std::string::npos) {
      amrex::AllPrint() << "FAIL: " << error.what() << '\n';
      status = 1;
    }
  }
  amrex::ParallelDescriptor::ReduceIntSum(caught);
  if (!expected.empty() && caught == 0) {
    amrex::Print() << "FAIL: expected rejection: " << expected << '\n';
    status = 1;
  }
  amrex::ParallelDescriptor::ReduceIntMax(status);
  amrex::Finalize();
  return status;
}
