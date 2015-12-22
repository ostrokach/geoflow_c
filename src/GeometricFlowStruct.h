#ifndef __GEOMETRICFLOWSTRUCT_H
#define __GEOMETRICFLOWSTRUCT_H

// Boundary element type
// TODO: "focus" (higher priority) and "sdh" should be added;  Only type
// implemented is MDH
enum BoundaryType{ ZERO, SDH, MDH, FOCUS, MAP};

//
//  input
//
struct GeometricFlowInput {

   enum BoundaryType m_boundaryCondition; // TODO: make this an enum?
   int m_vdwdispersion;  // 1/0, on/off 
   double m_gamma;
   double m_grid;
   double m_etolSolvation;
   double m_tol;
   double m_pdie;
   double m_sdie;
   double m_press;
   double m_bconc;

} ;

//
//  output
//
struct GeometricFlowOutput{

	double m_area,
		m_volume,
		m_attint,
		m_sumpot,
		m_totalSolvation,
		m_nonpolarSolvation,
		m_elecSolvation;

} ;

#endif
