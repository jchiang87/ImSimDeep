#include <fstream>
#include <sstream>

#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/geom/Point.h"
#include "desc/imsimdeep/InstcatTools.h"

namespace desc {
namespace imsimdeep {

   double InstcatTools::
   ang_sep(double ra0, double dec0, double ra1, double dec1) {
      lsst::afw::geom::Point2D p0(ra0, dec0);
      lsst::afw::geom::Point2D p1(ra1, dec1);
      lsst::afw::coord::Coord coord0(p0);
      lsst::afw::coord::Coord coord1(p1);
      return coord0.angularSeparation(coord1).asDegrees();
   }

   void InstcatTools::sky_cone_select(const std::string & infile,
                                      double ra, double dec, double radius,
                                      const std::string & outfile) {
      std::ifstream input(infile.c_str());
      std::ofstream output(outfile.c_str());
      std::string line;
      while (std::getline(input, line, '\n')) {
         if (line.substr(0, 6) != "object") {
            output << line << std::endl;
         } else {
            std::string command;
            std::string objectID;
            double ra_obj, dec_obj;
            std::istringstream ss;
            ss.str(line);
            ss >> command >> objectID >> ra_obj >> dec_obj;
            if (ang_sep(ra, dec, ra_obj, dec_obj) <= radius) {
               output << line << std::endl;
            }
         }
      }
      input.close();
      output.close();
   }
} // namespace imsimdeep
} // namespace desc
