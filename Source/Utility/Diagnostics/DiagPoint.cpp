#include "DiagPoint.H"
#include <AMReX_GpuLaunch.H>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

void
DiagPoint::init(const std::string& a_prefix, std::string_view a_diagName)
{
  amrex::ParmParse pp(a_prefix);
  int interval = -1;
  pp.query("int", interval);
  if (interval <= 0) {
    amrex::Abort("DiagPoint: a positive int is required for " + a_prefix);
  }
  amrex::Real period = -1.0;
  pp.query("per", period);
  if (period > 0.0) {
    amrex::Abort("DiagPoint: per is not supported for " + a_prefix);
  }
  if (pp.countval("filters") > 0) {
    amrex::Abort("DiagPoint: filters are not supported for " + a_prefix);
  }
  DiagBase::init(a_prefix, a_diagName);

  if (pp.countval("position") != AMREX_SPACEDIM) {
    amrex::Abort(
      "DiagPoint: position must have AMREX_SPACEDIM coordinates for " +
      a_prefix);
  }
  amrex::Vector<amrex::Real> position;
  pp.getarr("position", position, 0, AMREX_SPACEDIM);
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (!std::isfinite(position[d])) {
      amrex::Abort("DiagPoint: position must be finite for " + a_prefix);
    }
    m_position[d] = position[d];
  }
  if (pp.countval("field_names") == 0) {
    amrex::Abort("DiagPoint: field_names must not be empty for " + a_prefix);
  }
  pp.getarr("field_names", m_fieldNames);
  pp.query("volume_fraction_field", m_volumeFractionField);
}

void
DiagPoint::addVars(amrex::Vector<std::string>& a_varList)
{
  a_varList.insert(a_varList.end(), m_fieldNames.begin(), m_fieldNames.end());
  if (!m_volumeFractionField.empty()) {
    a_varList.push_back(m_volumeFractionField);
  }
}

void
DiagPoint::prepare(
  int a_nlevels,
  const amrex::Vector<amrex::Geometry>& a_geoms,
  const amrex::Vector<amrex::BoxArray>& a_grids,
  const amrex::Vector<amrex::DistributionMapping>& a_dmap,
  const amrex::Vector<std::string>& a_varNames)
{
  if (first_time) {
    amrex::Vector<int> indices;
    for (const auto& field : m_fieldNames) {
      indices.push_back(getFieldIndex(field, a_varNames));
    }
    if (!m_volumeFractionField.empty()) {
      indices.push_back(getFieldIndex(m_volumeFractionField, a_varNames));
    }
    m_fieldIndices_d.resize(indices.size());
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, indices.begin(), indices.end(),
      m_fieldIndices_d.begin());
    first_time = false;
  }

  AMREX_ALWAYS_ASSERT(a_nlevels > 0);
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (
      m_position[d] < a_geoms[0].ProbLo(d) ||
      m_position[d] > a_geoms[0].ProbHi(d)) {
      amrex::Abort(
        "DiagPoint: position is outside the physical domain for " + m_diagfile);
    }
  }
  // Recompute both level and owner: grids and their MPI mapping can change.
  m_level = -1;
  for (int lev = a_nlevels - 1; lev >= 0; --lev) {
    const auto& geom = a_geoms[lev];
    const auto& domain = geom.Domain();
    amrex::IntVect cell;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      const int offset = static_cast<int>(
        std::floor((m_position[d] - geom.ProbLo(d)) / geom.CellSize(d)));
      // Include the physical upper boundary in the last interior cell.
      cell[d] = std::clamp(
        domain.smallEnd(d) + offset, domain.smallEnd(d), domain.bigEnd(d));
    }
    const auto hits = a_grids[lev].intersections(amrex::Box(cell, cell));
    if (!hits.empty()) {
      m_level = lev;
      m_gridIndex = hits[0].first;
      m_owner = a_dmap[lev][m_gridIndex];
      m_cell = cell;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        m_cellCenter[d] =
          geom.ProbLo(d) +
          (cell[d] - domain.smallEnd(d) + 0.5) * geom.CellSize(d);
      }
      break;
    }
  }
  if (m_level < 0) {
    amrex::Abort("DiagPoint: no valid cell covers position for " + m_diagfile);
  }
}

void
DiagPoint::processDiag(
  int a_nstep,
  const amrex::Real& a_time,
  const amrex::Vector<const amrex::MultiFab*>& a_state,
  const amrex::Vector<std::string>& /*a_varNames*/)
{
  AMREX_ALWAYS_ASSERT(m_level >= 0 && m_level < a_state.size());
  const int count = static_cast<int>(m_fieldIndices_d.size());
  amrex::Vector<amrex::Real> values(count, 0.0);
  if (amrex::ParallelDescriptor::MyProc() == m_owner) {
    const auto data = a_state[m_level]->const_array(m_gridIndex);
    const auto cell = m_cell;
    const auto* indices = m_fieldIndices_d.data();
    amrex::Gpu::DeviceVector<amrex::Real> sampled(count);
    auto* sample = sampled.data();
    amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE(int n) noexcept {
      sample[n] = data(cell, indices[n]);
    });
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, sampled.begin(), sampled.end(), values.begin());
  }
  amrex::ParallelDescriptor::ReduceRealSum(
    values.data(), count, amrex::ParallelDescriptor::IOProcessorNumber());
  if (amrex::ParallelDescriptor::IOProcessor()) {
    if (!m_volumeFractionField.empty() && !(values.back() > 0.0)) {
      amrex::Abort(
        "DiagPoint: selected cell is fully covered by EB for " + m_diagfile);
    }
    writeData(a_nstep, a_time, values);
  }
}

std::string
DiagPoint::positionHeader() const
{
  std::ostringstream header;
  header << "# DiagPoint position" << std::scientific
         << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10);
  for (const auto coordinate : m_position) {
    header << ' ' << coordinate;
  }
  return header.str();
}

std::string
DiagPoint::columnsHeader() const
{
  const char* axes[] = {"x", "y", "z"};
  std::ostringstream header;
  header << "# step time";
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    header << ' ' << axes[d];
  }
  header << " level";
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    header << " cell_" << axes[d];
  }
  for (const auto& field : m_fieldNames) {
    header << ' ' << field;
  }
  return header.str();
}

void
DiagPoint::initOutput()
{
  std::error_code error;
  const bool exists = std::filesystem::exists(m_diagfile, error);
  if (error) {
    amrex::Abort("DiagPoint: cannot open " + m_diagfile);
  }
  bool needsHeader = true;
  if (exists) {
    std::ifstream input(m_diagfile, std::ios::binary);
    if (!input) {
      amrex::Abort("DiagPoint: cannot open " + m_diagfile);
    }
    if (input.peek() != std::ifstream::traits_type::eof()) {
      needsHeader = false;
      input.seekg(-1, std::ios::end);
      char last = '\0';
      input.get(last);
      if (last != '\n') {
        amrex::Abort("DiagPoint: invalid history record in " + m_diagfile);
      }
      input.seekg(0);
      std::string line;
      std::getline(input, line);
      if (line != positionHeader()) {
        amrex::Abort(
          "DiagPoint: existing position does not match " + m_diagfile);
      }
      std::getline(input, line);
      if (line != columnsHeader()) {
        amrex::Abort("DiagPoint: existing columns do not match " + m_diagfile);
      }
      const int remaining =
        1 + 2 * AMREX_SPACEDIM + static_cast<int>(m_fieldNames.size());
      while (std::getline(input, line)) {
        std::istringstream record(line);
        int step = -1;
        amrex::Real time = 0.0;
        bool valid = static_cast<bool>(record >> step >> time);
        for (int col = 0; col < remaining; ++col) {
          amrex::Real value;
          valid = valid && static_cast<bool>(record >> value);
        }
        std::string extra;
        valid = valid && !(record >> extra) && step >= 0 && std::isfinite(time);
        if (
          !valid ||
          (m_hasRecords && (step <= m_lastStep || time <= m_lastTime))) {
          amrex::Abort("DiagPoint: invalid history record in " + m_diagfile);
        }
        m_lastStep = step;
        m_lastTime = time;
        m_hasRecords = true;
      }
      if (!input.eof()) {
        amrex::Abort("DiagPoint: cannot read " + m_diagfile);
      }
    }
  }
  if (needsHeader) {
    std::ofstream output(m_diagfile, std::ios::app | std::ios::binary);
    output << positionHeader() << '\n' << columnsHeader() << '\n';
    output.flush();
    if (!output) {
      amrex::Abort("DiagPoint: cannot open or write " + m_diagfile);
    }
  }
  m_outputInitialized = true;
}

void
DiagPoint::writeData(
  int a_nstep, amrex::Real a_time, const amrex::Vector<amrex::Real>& a_values)
{
  if (!m_outputInitialized) {
    initOutput();
  }
  if (m_hasRecords) {
    if (a_nstep == m_lastStep && a_time == m_lastTime) {
      return;
    }
    if (a_nstep < m_lastStep || a_time < m_lastTime) {
      amrex::Abort(
        "DiagPoint: output is older than the last record in " + m_diagfile +
        "; use a new file when restarting from an earlier checkpoint");
    }
    if (a_nstep == m_lastStep || a_time == m_lastTime) {
      amrex::Abort("DiagPoint: inconsistent step/time for " + m_diagfile);
    }
  }
  std::ofstream output(m_diagfile, std::ios::app | std::ios::binary);
  output << std::scientific
         << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10)
         << a_nstep << ' ' << a_time;
  for (const auto coordinate : m_position) {
    output << ' ' << coordinate;
  }
  output << ' ' << m_level;
  for (const auto coordinate : m_cellCenter) {
    output << ' ' << coordinate;
  }
  for (int n = 0; n < m_fieldNames.size(); ++n) {
    output << ' ' << a_values[n];
  }
  output << '\n';
  output.flush();
  if (!output) {
    amrex::Abort("DiagPoint: cannot open or write " + m_diagfile);
  }
  m_lastStep = a_nstep;
  m_lastTime = a_time;
  m_hasRecords = true;
}
