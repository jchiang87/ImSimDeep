#ifndef desc_imsimdeep_InstcatTools_h
#define desc_imsimdeep_InstcatTools_h

#include <string>

namespace desc {
   namespace imsimdeep {

      class InstcatTools {
      public:
#ifndef SWIG
         InstcatTools() {}
#endif
         static double ang_sep(double ra0, double dec0,
                               double ra1, double dec1);

         static void sky_cone_select(const std::string & infile,
                                     double ra, double dec, double radius,
                                     const std::string & outfile);
      private:
      };

   } // namespace imsimdeep
} // namespace desc

#endif // desc_imsimdeep_InstcatTools_h
