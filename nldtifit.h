// Nonlinear estimation of diffusion tensor
// using a Gaussian or Rician noise model
//
// nldtifit 
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2011 University of Oxford 
//

namespace NLDTIFIT {

class nldtifitException: public std::exception
{
private:
  std::string m_msg;
public:
  nldtifitException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("nldtifit: msg=" + m_msg).c_str();
  }

  ~nldtifitException() throw() {}
};

enum NbHood {NONE = 1, FACE=7, EDGE=19, CORNER=27};
enum NoiseModel {GAUSSIAN=1, RICIAN=2};

NEWIMAGE::volume4D<float> nldtifit(// Data, b=0 and diffusion weighted
                                   const NEWIMAGE::volume4D<float>&   scans,
                                   // Initial guesses of the parameters in the order s2, s0, R1, R2, R3, L1, L2, L3.
                                   // N.B. that s2, s0 and L* must all be positive
                                   const NEWIMAGE::volume4D<float>&   inpar,
                                   // Mask indicating where we want calculations to be performed
                                   const NEWIMAGE::volume<char>&      mask,
                                   // nx3 matrix with gradient directions in same order as in scans
                                   const NEWMAT::Matrix&              bvec,
                                   // nx1 vector with b-values in same order as in scans. N.B. that wa can have many different b-values
                                   const NEWMAT::ColumnVector&        bval,
                                   // Indicates what neighbourhood of voxels we want to base s2 estimate on.
                                   NbHood                             nbhood,
                                   // Gaussian or Rician
                                   NoiseModel                         nm);

void verify_input(const NEWIMAGE::volume4D<float>&   scans,
                  const NEWIMAGE::volume4D<float>&   inpar,
                  const NEWIMAGE::volume<char>&      mask,
                  const NEWMAT::Matrix&              bvec,
                  const NEWMAT::ColumnVector&        bval,
                  NbHood                             nbhood,
                  NoiseModel                         nm);

void deal_results(// Input
                  const double                      *theta,
                  int                               ii,
                  int                               jj,
                  int                               kk,
                  // Output
                  NEWIMAGE::volume4D<float>&        ovol);

void repack_g_and_b(// Input 
                    const NEWMAT::Matrix&         bvec,
                    const NEWMAT::ColumnVector&   bval,
		    // Output
		    double                        *g,
                    double                        *b);

int load_data_and_initial_guesses(// Input
                                  const NEWIMAGE::volume4D<float>&  scans,
                                  const NEWIMAGE::volume4D<float>&  inpar,
                                  const NEWIMAGE::volume<char>&     mask,
                                  int                               ii,
                                  int                               jj,
                                  int                               kk,
                                  NbHood                            nbhood,
                                  // Output
                                  double                            *y,
                                  double                            *ig);

double *load_data_single_voxel(const NEWIMAGE::volume4D<float>&  scans,
                               int                               ii,
                               int                               jj,
                               int                               kk,
                               // Output
                               double                            *y);

double *load_guess_single_voxel(const NEWIMAGE::volume4D<float>&  inpar,
                                int                               ii,
                                int                               jj,
                                int                               kk,
                                // Output
                                double                            *ig);

int fit_RLR(/* Input */
            double   *g,      /* Gradients, nx3 Matlab matrix. */
            double   *b,      /* b-values. */
            double   *y,      /* Data. */
            int      n,       /* # of data points. */
            int      nvox,    /* # of voxels. */
            double   *pmu,    /* Means of priors. */
            double   *pvar,   /* Variances of priors. */
            int      *pi,     /* Indicies (into theta) of priors. */
            int      np,      /* # of parameters with priors. */
            int      nm,      /* Noise model 1->Gauss, 2->Rician. */
            /* Input/Output */
            double   *theta,  /* Parameters. */
            double   *ei,     /* Expected intensities. */
            double   *hist);  /* History of log-likelihoods. */

} // end namespace NLDTIFIT
