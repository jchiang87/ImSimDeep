#ifndef desc_imsimdeep_instcat_utils_h
#define desc_imsimdeep_instcat_utils_h

#include <string>

namespace desc {
   namespace imsimdeep {

      double ang_sep(double ra0, double dec0, double ra1, double dec1);

      void sky_cone_select(const std::string & infile,
                           double ra, double dec, double radius,
                           const std::string & outfile);
   } // namespace imsimdeep
} // namespace desc

#endif // desc_imsimdeep_instcat_utils_h
