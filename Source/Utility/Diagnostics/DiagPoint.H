#ifndef DIAGPOINT_H
#define DIAGPOINT_H

#include "DiagBase.H"
#include <AMReX_Geometry.H>

class DiagPoint : public DiagBase::Register<DiagPoint>
{
public:
  static std::string identifier() { return "DiagPoint"; }

  void init(const std::string& a_prefix, std::string_view a_diagName) override;
  void addVars(amrex::Vector<std::string>& a_varList) override;
  void prepare(
    int a_nlevels,
    const amrex::Vector<amrex::Geometry>& a_geoms,
    const amrex::Vector<amrex::BoxArray>& a_grids,
    const amrex::Vector<amrex::DistributionMapping>& a_dmap,
    const amrex::Vector<std::string>& a_varNames) override;
  void processDiag(
    int a_nstep,
    const amrex::Real& a_time,
    const amrex::Vector<const amrex::MultiFab*>& a_state,
    const amrex::Vector<std::string>& a_varNames) override;
  void close() override {}

private:
  std::string positionHeader() const;
  std::string columnsHeader() const;
  void initOutput();
  void writeData(
    int a_nstep,
    amrex::Real a_time,
    const amrex::Vector<amrex::Real>& a_values);

  amrex::RealArray m_position{};
  amrex::RealArray m_cellCenter{};
  amrex::IntVect m_cell{};
  int m_level{-1};
  int m_gridIndex{-1};
  int m_owner{-1};
  amrex::Vector<std::string> m_fieldNames;
  std::string m_volumeFractionField;
  // If requested, the final component rejects fully covered EB cells.
  amrex::Gpu::DeviceVector<int> m_fieldIndices_d;
  bool m_outputInitialized{false};
  bool m_hasRecords{false};
  int m_lastStep{-1};
  amrex::Real m_lastTime{0.0};
};

#endif
