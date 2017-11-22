//Common items for the Kepler code

using System;
namespace ASCOM.Astrometry.Kepler {
   sealed class KeplerGlobalCode
   {

      #region Private Structures
      internal struct plantbl
      {
         internal int maxargs;
         internal int[] max_harmonic;
         internal int max_power_of_t;
         internal int[] arg_tbl;
         internal double[] lon_tbl;
         internal double[] lat_tbl;
         internal double[] rad_tbl;
         internal double distance;
         internal double timescale;
         internal double trunclvl;

         internal plantbl(int ma, int[] mh, int mpt, int[] at, double[] lot, double[] lat, double[] rat, double dis, double ts, double tl)
         {
            maxargs = ma;
            max_harmonic = mh;
            max_power_of_t = mpt;
            arg_tbl = at;
            lon_tbl = lot;
            lat_tbl = lat;
            rad_tbl = rat;
            distance = dis;
            timescale = ts;
            trunclvl = tl;
         }
      }

      internal struct orbit
      {
         internal string obname; ///* name of the object */
         internal double epoch; ///* epoch of orbital elements */
         internal double i; ///* inclination	*/
         internal double W; ///* longitude of the ascending node */
         internal double wp; ///* argument of the perihelion */
         internal double a; ///* mean distance (semimajor axis) */
         internal double dm; ///* daily motion */
         internal double ecc; ///* eccentricity */
         internal double M; ///* mean anomaly */
         internal double equinox; ///* epoch of equinox and ecliptic */
         internal double mag; ///* visual magnitude at 1AU from earth and sun */
         internal double sdiam; ///* equatorial semidiameter at 1au, arc seconds */
                                ///* The following used by perterbation formulas: */
         internal plantbl ptable;
         internal double L; ///* computed mean longitude */
         internal double r; ///* computed radius vector */
         internal double plat; ///* perturbation in ecliptic latitude */

         internal orbit(string obn, double ep, double i_p, double W_p, double wp_p,
            double a_p, double dm_p, double ecc_p, double M_p,
            double eq, double mg, double sd, plantbl pt, double L_p,
            double r_p, double pl)
         {
            obname = obn;
            epoch = ep;
            i = i_p;
            W = W_p;
            wp = wp_p;
            a = a_p;
            dm = dm_p;
            ecc = ecc_p;
            M = M_p;
            equinox = eq;
            mag = mg;
            sdiam = sd;
            ptable = pt;
            L = L_p;
            r = r_p;
            plat = pl;
         }
      }
      #endregion

      #region Constants
      internal const int NARGS = 18;

      ///* Conversion factors between degrees and radians */
      private const double DTR = 0.017453292519943295;
      private const double RTD = 57.295779513082323;
      private const double RTS = 206264.80624709636; ///* arc seconds per radian */
      private const double STR = 0.00000484813681109536; ///* radians per arc second */
      private const double PI = 3.1415926535897931;
      private const double TPI = 2.0 * PI;

      ///* Standard epochs.  Note Julian epochs (J) are measured in
      // * years of 365.25 days.
      // */
      private const double J2000 = 2451545.0; ///* 2000 January 1.5 */
      private const double B1950 = 2433282.423; ///* 1950 January 0.923 Besselian epoch */
      private const double J1900 = 2415020.0; ///* 1900 January 0, 12h UT */

      ///* Constants used elsewhere. These are DE403 values. */
      private const double aearth = 6378137.0; ///* Radius of the earth, in meters.  */
      private const double au = 149597870.691; ///* Astronomical unit, in kilometers.  */
      private const double emrat = 81.300585; ///* Earth/Moon mass ratio.  */
      private const double Clight = 299792.458; ///* Speed of light, km/sec  */
      // private const double Clightaud = null; ///* C in au/day  */
      #endregion

      #region Utility Routines
      //// ----------------
      //// Utility routines
      //// ----------------

      //// Obliquity of the ecliptic at Julian date J
      //// according to the DE403 values. Refer to
      //// S. Moshier's aa54e sources.

      internal static void epsiln(double J, ref double eps, ref double coseps, ref double sineps)
      {
         double T = 0;

         T = System.Convert.ToDouble((J - 2451545.0) / 365250.0); //// T / 10
         eps = System.Convert.ToDouble(((((((((((0.000000000245 * T + 0.00000000579) * T + 0.0000002787) * T
            + 0.000000712) * T - 0.00003905) * T - 0.0024967) * T
            - 0.005138) * T + 1.9989) * T - 0.0175) * T - 468.3396) * T
            + 84381.406173) * STR);
         coseps = System.Convert.ToDouble(Math.Cos(eps));
         sineps = System.Convert.ToDouble(Math.Sin(eps));
      }

      ///* Precession of the equinox and ecliptic
      // * from epoch Julian date J to or from J2000.0
      // *
      // * Program by Steve Moshier.  */
      //
      ///* James G. Williams, "Contributions to the Earth's obliquity rate,
      //   precession, and nutation,"  Astron. J. 108, 711-724 (1994)  */

      // /* Corrections to Williams (1994) introduced in DE403.  */
      internal static double[] pAcof = new double[] {-0.000000000866, -0.00000004759, 0.0000002424, 0.000013095,
      0.00017451, -0.0018055, -0.235316, 0.076, 110.5414, 50287.91959};

      internal static double[] nodecof = new double[] {0.00000000000000066402, -0.00000000000000269151, -0.000000000001547021,
      0.000000000007521313, 0.00000000019, -0.00000000354, -0.00000018103,
      0.000000126, 0.00007436169, -0.04207794833, 3.052115282424};

      internal static double[] inclcof = new double[] {0.00000000000000012147, System.Convert.ToDouble(7.3759E-17), -0.0000000000000826287,
      0.000000000000250341, 0.000000000024650839, -0.000000000054000441,
      0.00000000132115526, -0.0000006012, -0.0000162442, 0.00227850649, 0.0};

      ///* Subroutine arguments:
      // *
      // * R = rectangular equatorial coordinate vector to be precessed.
      // *     The result is written back into the input vector.
      // * J = Julian date
      // * direction =
      // *      Precess from J to J2000: direction = 1
      // *      Precess from J2000 to J: direction = -1
      // * Note that if you want to precess from J1 to J2, you would
      // * first go from J1 to J2000, then call the program again
      // * to go from J2000 to J2.
      // */

      internal static void precess(ref double[] R, double J, int direction)
      {
         double A = 0;
         double B = 0;
         double T = 0;
         double pA = 0;
         double W = 0;
         double z = 0;
         double[] x = new double[4];
         double[] p = null;
         double eps = 0;
         double coseps = 0;
         double sineps = 0;
         int i = 0;

         if (J == J2000) {
            return;
         }
         ///* Each precession angle is specified by a polynomial in
         // * T = Julian centuries from J2000.0.  See AA page B18.
         // */
         T = System.Convert.ToDouble((J - J2000) / 36525.0);

         ///* Implementation by elementary rotations using Laskar's expansions.
         // * First rotate about the x axis from the initial equator
         // * to the ecliptic. (The input is equatorial.)
         // */
         if (direction == 1) {
            epsiln(J, ref eps, ref coseps, ref sineps); ///* To J2000 */
         }
         else {
            epsiln(J2000, ref eps, ref coseps, ref sineps); ///* From J2000 */
         }
         x[0] = R[0];
         z = coseps * R[1] + sineps * R[2];
         x[2] = System.Convert.ToDouble(-sineps * R[1] + coseps * R[2]);
         x[1] = z;

         ///* Precession in longitude	 */
         T /= 10.0; ///* thousands of years */
         p = pAcof;
         pA = p[0];
         for (i = 1; i <= 9; i++) {
            pA = pA * T + p[i];
         }
         pA *= STR * T;

         ///* Node of the moving ecliptic on the J2000 ecliptic.*/
         p = nodecof;
         W = p[0];
         for (i = 1; i <= 10; i++) {
            W = W * T + p[i];
         }
         ///* Rotate about z axis to the node.*/
         if (direction == 1) {
            z = W + pA;
         }
         else {
            z = W;
         }
         B = System.Convert.ToDouble(Math.Cos(z));
         A = System.Convert.ToDouble(Math.Sin(z));
         z = B * x[0] + A * x[1];
         x[1] = System.Convert.ToDouble(-A * x[0] + B * x[1]);
         x[0] = z;

         ///* Rotate about new x axis by the inclination of the moving
         // * ecliptic on the J2000 ecliptic.
         // */
         p = inclcof;
         z = p[0];
         for (i = 1; i <= 10; i++) {
            z = z * T + p[i];
         }
         if (direction == 1) {
            z = System.Convert.ToDouble(-z);
         }
         B = System.Convert.ToDouble(Math.Cos(z));
         A = System.Convert.ToDouble(Math.Sin(z));
         z = B * x[1] + A * x[2];
         x[2] = System.Convert.ToDouble(-A * x[1] + B * x[2]);
         x[1] = z;

         ///* Rotate about new z axis back from the node.	 */
         if (direction == 1) {
            z = System.Convert.ToDouble(-W);
         }
         else {
            z = System.Convert.ToDouble(-W - pA);
         }
         B = System.Convert.ToDouble(Math.Cos(z));
         A = System.Convert.ToDouble(Math.Sin(z));
         z = B * x[0] + A * x[1];
         x[1] = System.Convert.ToDouble(-A * x[0] + B * x[1]);
         x[0] = z;

         ///* Rotate about x axis to final equator.	 */
         if (direction == 1) {
            epsiln(J2000, ref eps, ref coseps, ref sineps);
         }
         else {
            epsiln(J, ref eps, ref coseps, ref sineps);
         }
         z = coseps * x[1] - sineps * x[2];
         x[2] = sineps * x[1] + coseps * x[2];
         x[1] = z;

         for (i = 0; i <= 2; i++) {
            R[i] = x[i];
         }
      }

      internal static double atan4(double x, double y)
      {

         double z = 0;
         double w = 0;
         int code = 0;

         code = 0;

         if (x < 0.0) {
            code = 2;
         }
         if (y < 0.0) {
            code = code | 1;
         }

         if (x == 0.0) {
            if ((code & 1) > 0) {
               return (1.5 * PI);
            }
            if (y == 0.0) {
               return (0.0);
            }
            return (0.5 * PI);
         }

         if (y == 0.0) {
            if ((code & 2) > 0) {
               return (PI);
            }
            return (0.0);
         }

         switch (code) {
            case 0:
               w = 0.0;
               break;
            case 1:
               w = System.Convert.ToDouble(2.0 * PI);
               break;
            case 2:
               break;
            case 3:
               w = PI;
               break;
            default:
               break;
         }

         z = System.Convert.ToDouble(Math.Atan(y / x));

         return (w + z);
      }

      ////
      //// Reduce x modulo 2 pi
      ////
      internal static double modtp(double x)
      {

         double y = 0;

         y = System.Convert.ToDouble(Math.Floor(x / TPI));
         y = x - y * TPI;
         while (y < 0.0) {
            y += TPI;
         }
         while (y >= TPI) {
            y -= TPI;
         }
         return (y);
      }

      ////
      ////  Reduce x modulo 360 degrees
      ////
      internal static double mod360(double x)
      {

         int k = 0;
         double y = 0;

         k = (int)(x / 360.0);
         y = x - k * 360.0;
         while (y < 0.0) {
            y += 360.0;
         }
         while (y > 360.0) {
            y -= 360.0;
         }
         return (y);
      }


      ///* Program to solve Keplerian orbit
      // * given orbital parameters and the time.
      // * Returns Heliocentric equatorial rectangular coordinates of
      // * the object.
      // *
      // * This program detects several cases of given orbital elements.
      // *
      // * If a program for perturbations is pointed to, it is called
      // * to calculate all the elements.
      // *
      // * If there is no program, then the mean longitude is calculated
      // * from the mean anomaly and daily motion.
      // *
      // * If the daily motion is not given, it is calculated
      // * by Kepler's law.
      // *
      // * If the eccentricity is given to be 1.0, it means that
      // * meandistance is really the perihelion distance, as in a comet
      // * specification, and the orbit is parabolic.
      // *
      // * Reference: Taff, L.G., "Celestial Mechanics, A Computational
      // * Guide for the Practitioner."  Wiley, 1985.
      // */


      internal static void KeplerCalc(double J, orbit e, ref double[] rect)
      {

         double[] polar = new double[4];
         double alat = 0;
         double E1 = 0;
         double M = 0;
         double W = 0;
         double temp = 0;
         double epoch = 0;
         double inclination = 0;
         double ascnode = 0;
         double argperih = 0;
         double meandistance = 0;
         double dailymotion = 0;
         double eccent = 0;
         double meananomaly = 0;
         double r = 0;
         double coso = 0;
         double sino = 0;
         double cosa = 0;
         double eps = 0;
         double coseps = 0;
         double sineps = 0;
         //Dim TL As New TraceLogger("", "KeplerCalc")
         //TL.Enabled = True
         //TL.LogMessage("KepCalc", "Started")
         ////
         //// Call program to compute position, if one is supplied.
         ////
         if (e.ptable.lon_tbl[0] != 0.0) {
            if (e.obname == "Earth") {
               //TL.LogMessage("KepCalc", "Before G3Plan Earth")
               g3plan(J, e.ptable, ref polar, 3);
            }
            else {
               //TL.LogMessage("KepCalc", "Before G3Plan Not Earth")
               gplan(J, e.ptable, ref polar);
            }
            E1 = polar[0]; ///* longitude */
            e.L = E1;
            W = polar[1]; ///* latitude */
            r = polar[2]; ///* radius */
            e.r = r;
            e.epoch = J;
            e.equinox = J2000;
            //TL.LogMessage("KepCalc", "After G3Plan")
            goto kepdon;
         }

         //// -----------------------------
         //// Compute from orbital elements
         //// -----------------------------
         //TL.LogMessage("KepCalc", "Compute orbital elements")
         e.equinox = J2000; //// Always J2000 coordinates
         epoch = e.epoch;
         inclination = e.i;
         ascnode = e.W * DTR;
         argperih = e.W;
         meandistance = e.a; ///* semimajor axis */
         dailymotion = e.dm;
         eccent = e.ecc;
         meananomaly = e.M;

         //// ---------
         //// Parabolic
         //// ---------
         if (eccent == 1.0) {
            //TL.LogMessage("KepCalc", "eccent=1.0")
            ////
            //// meandistance = perihelion distance, q
            //// epoch = perihelion passage date
            ////
            temp = meandistance * Math.Sqrt(meandistance);
            W = System.Convert.ToDouble((J - epoch) * 0.0364911624 / temp);
            ////
            //// The constant above is 3 k / sqrt(2),
            //// k = Gaussian gravitational constant = 0.01720209895
            ////
            E1 = 0.0;
            M = 1.0;
            while (Math.Abs(M) > 0.00000000001) {
               temp = E1 * E1;
               temp = System.Convert.ToDouble((2.0 * E1 * temp + W) / (3.0 * (1.0 + temp)));
               M = temp - E1;
               if (temp != 0.0) {
                  M /= temp;
               }
               E1 = temp;
            }
            r = meandistance * (1.0 + E1 * E1);
            M = System.Convert.ToDouble(Math.Atan(E1));
            M = System.Convert.ToDouble(2.0 * M);
            alat = M + (DTR * argperih);
            //// ----------
            //// Hyperbolic
            //// ----------
         }
         else if (eccent > 1.0) {
            //TL.LogMessage("KepCalc", "eccent > 1.0")
            ////
            //// The equation of the hyperbola in polar coordinates r, theta
            //// is r = a(e^2 - 1)/(1 + e cos(theta)) so the perihelion
            //// distance q = a(e-1), the "mean distance"  a = q/(e-1).
            ////
            meandistance = meandistance / (eccent - 1.0);
            temp = meandistance * Math.Sqrt(meandistance);
            W = System.Convert.ToDouble((J - epoch) * 0.01720209895 / temp);
            ///* solve M = -E + e sinh E */
            E1 = W / (eccent - 1.0);
            M = 1.0;
            while (Math.Abs(M) > 0.00000000001) {

               M = System.Convert.ToDouble(-E1 + eccent * Math.Sinh(E1) - W);
               E1 += M / (1.0 - eccent * Math.Cosh(E1));
            }
            r = meandistance * (-1.0 + eccent * Math.Cosh(E1));
            temp = System.Convert.ToDouble((eccent + 1.0) / (eccent - 1.0));
            M = System.Convert.ToDouble(Math.Sqrt(temp) * Math.Tanh(0.5 * E1));
            M = System.Convert.ToDouble(2.0 * Math.Atan(M));
            alat = M + (DTR * argperih);

            //// -----------
            //// Ellipsoidal
            //// -----------
         }
         else //		// if(ecc < 1)
         {
            //TL.LogMessage("KepCalc", "Ellipsoidal")
            ////
            //// Calculate the daily motion, if it is not given.
            ////
            if (dailymotion == 0.0) {

               ////
               //// The constant is 180 k / pi, k = Gaussian gravitational
               //// constant. Assumes object in heliocentric orbit is
               //// massless.
               ////
               dailymotion = System.Convert.ToDouble(0.9856076686 / (e.a * Math.Sqrt(e.a)));
            }
            dailymotion *= J - epoch;
            ////
            //// M is proportional to the area swept out by the radius
            //// vector of a circular orbit during the time between
            //// perihelion passage and Julian date J.
            //// It is the mean anomaly at time J.
            ////
            M = DTR * (meananomaly + dailymotion);
            M = modtp(M);
            ////
            //// If mean longitude was calculated, adjust it also
            //// for motion since epoch of elements.
            ////
            if ((e.L) != 0.0) {
               e.L += dailymotion;
               e.L = mod360(e.L);
            }
            ////
            //// By Kepler's second law, M must be equal to
            //// the area swept out in the same time by an
            //// elliptical orbit of same total area.
            //// Integrate the ellipse expressed in polar coordinates
            ////     r = a(1-e^2)/(1 + e cosW)
            //// with respect to the angle W to get an expression for the
            //// area swept out by the radius vector.  The area is given
            //// by the mean anomaly; the angle is solved numerically.
            ////
            //// The answer is obtained in two steps.  We first solve
            //// Kepler's equation
            ////    M = E - eccent*sin(E)
            //// for the eccentric anomaly E.  Then there is a
            //// closed form solution for W in terms of E.
            ////
            E1 = M; ///* Initial guess is same as circular orbit. */
            temp = 1.0;
            do {
               //// The approximate area swept out in the ellipse
               //// ...minus the area swept out in the circle
               temp = E1 - eccent * Math.Sin(E1) - M;
               //// ...should be zero.  Use the derivative of the error
               ////to converge to solution by Newton's method.
               E1 -= temp / (1.0 - (eccent * Math.Cos(E1)));
            } while (Math.Abs(temp) > 0.00000000001);

            ////
            //// The exact formula for the area in the ellipse is
            ////    2.0*atan(c2*tan(0.5*W)) - c1*eccent*sin(W)/(1+e*cos(W))
            //// where
            ////    c1 = sqrt( 1.0 - eccent*eccent )
            ////    c2 = sqrt( (1.0-eccent)/(1.0+eccent) ).
            //// Substituting the following value of W
            //// yields the exact solution.
            ////
            temp = System.Convert.ToDouble(Math.Sqrt((1.0 + eccent) / (1.0 - eccent)));
            W = System.Convert.ToDouble(2.0 * Math.Atan(temp * Math.Tan(0.5 * E1)));

            ////
            //// The true anomaly.
            ////
            W = modtp(W);

            meananomaly *= DTR;
            ////
            //// Orbital longitude measured from node
            //// (argument of latitude)
            ////
            if (e.L != 0.0) //// Mean longitude given
            {
               alat = System.Convert.ToDouble(((e.L) * DTR) + W - meananomaly - ascnode);
            }
            else {
               alat = W + (DTR * argperih); //// Mean longitude not given
            }
            ////
            //// From the equation of the ellipse, get the
            //// radius from central focus to the object.
            ////
            r = meandistance * (1.0 - eccent * eccent) / (1.0 + eccent * Math.Cos(W));
         }
         //TL.LogMessage("KepCalc", "Before All orbits")
         inclination *= DTR; //// Convert inclination to radians

         //// ----------
         //// ALL ORBITS
         //// ----------
         ////
         //// At this point:
         ////
         ////		alat		= argument of latitude (rad)
         ////		inclination	= inclination (rad)
         ////		r			= radius from central focus
         ////
         //// The heliocentric ecliptic longitude of the objectis given by:
         ////
         ////   tan(longitude - ascnode)  =  cos(inclination) * tan(alat)
         ////
         coso = System.Convert.ToDouble(Math.Cos(alat));
         sino = System.Convert.ToDouble(Math.Sin(alat));
         W = sino * Math.Cos(inclination);
         E1 = System.Convert.ToDouble(atan4(coso, W) + ascnode);

         ////
         //// The ecliptic latitude of the object
         ////
         W = System.Convert.ToDouble(Math.Asin(sino * Math.Sin(inclination)));

         //// ------------------------------------
         //// Both from DE404 and from elements...
         //// ------------------------------------
         ////
         //// At this point we have the heliocentric ecliptic polar
         //// coordinates of the body.
         ////
         kepdon:

         ////
         //// Convert to heliocentric ecliptic rectangular coordinates,
         //// using the perturbed latitude.
         ////
         rect[2] = r * Math.Sin(W);
         cosa = System.Convert.ToDouble(Math.Cos(W));
         rect[1] = r * cosa * Math.Sin(E1);
         rect[0] = r * cosa * Math.Cos(E1);

         ////
         //// Convert from heliocentric ecliptic rectangular
         //// to heliocentric equatorial rectangular coordinates
         //// by rotating epsilon radians about the x axis.
         ////
         //TL.LogMessage("KepCalc", "Before epsiln")
         epsiln(e.equinox, ref eps, ref coseps, ref sineps);
         W = coseps * rect[1] - sineps * rect[2];
         M = sineps * rect[1] + coseps * rect[2];
         rect[1] = W;
         rect[2] = M;

         ////
         //// Precess the equatorial (rectangular) coordinates to the
         //// ecliptic & equinox of J2000.0, if not already there.
         ////
         //TL.LogMessage("KepCalc", "Before precess")
         precess(ref rect, e.equinox, 1);

         ////
         //// If earth, adjust from earth-moon barycenter to earth
         //// by AA page E2.
         ////
         //TL.LogMessage("KepCalc", "Before embofs")
         if (e.obname == "Earth") {
            embofs(J, ref rect, ref r); ///* see embofs() below */
         }
         //TL.LogMessage("KepCalc", "Exited")

      }

      ////
      //// Adjust position from Earth-Moon barycenter to Earth
      ////
      //// J = Julian day number
      //// emb = Equatorial rectangular coordinates of EMB.
      //// pr = Earth's distance to the Sun (au)
      ////
      internal static void embofs(double J, ref double[] ea, ref double pr)
      {

         double[] pm = new double[4];
         double[] polm = new double[4];
         double a = 0;
         double b = 0;
         int i = 0;

         //Dim TL As New TraceLogger("", "Embofs")
         //TL.Enabled = True
         //TL.LogMessage("Embofs", "Start")
         ////
         //// Compute the vector Moon - Earth.
         ////
         //TL.LogMessage("Embofs", "Before GMoon")
         gmoon(J, ref pm, ref polm);
         //TL.LogMessage("Embofs", "After GMoon")

         ////
         //// Precess the lunar position
         //// to ecliptic and equinox of J2000.0
         ////
         precess(ref pm, J, 1);
         //TL.LogMessage("Embofs", "After Precess")

         ////
         //// Adjust the coordinates of the Earth
         ////
         a = System.Convert.ToDouble(1.0 / (emrat + 1.0));
         b = 0.0;
         for (i = 0; i <= 2; i++) {
            ea[i] = System.Convert.ToDouble(ea[i] - a * pm[i]);
            b = b + ea[i] * ea[i];
         }

         ////
         //// Sun-Earth distance.
         ////
         pr = System.Convert.ToDouble(Math.Sqrt(b));
      }
      #endregion

      #region MajElems
      ///* Orbits for each planet.  The indicated orbital elements are
      //* not actually used, since the positions are are now calculated
      //* from a formula.  Magnitude and semidiameter are still used.
      //*/


      ///* January 5.0, 1987 */
      internal static orbit mercury = new orbit("Mercury", 2446800.5, 7.0048, 48.177, 29.074, 0.387098, 4.09236,
         0.205628, 198.7199, 2446800.5, -0.42, 3.36, mer404, 0.0, 0.0, 0.0);

      ///* Note the calculated apparent visual magnitude for Venus is not very accurate. */
      internal static orbit venus = new orbit("Venus", 2446800.5, 3.3946, 76.561, 54.889, 0.723329, 1.60214, 0.006757, 9.0369, 2446800.5, -4.4, 8.34, ven404, 0.0, 0.0, 0.0);

      ///* Fixed numerical values will be used for earth if read in from a file named earth.orb.  See kfiles.c, kep.h. */
      internal static orbit earthplanet = new orbit("Earth", 2446800.5, 0.0, 0.0, 102.884, 0.999999, 0.985611, 0.016713, 1.1791, 2446800.5, -3.86, 0.0, ear404, 0.0, 0.0, 0.0);

      internal static orbit mars = new orbit("Mars", 2446800.5, 1.8498, 49.457, 286.343, 1.52371, 0.524023, 0.093472, 53.1893, 2446800.5, -1.52, 4.68, mar404, 0.0, 0.0, 0.0);

      internal static orbit jupiter = new orbit("Jupiter", 2446800.5, 1.3051, 100.358, 275.129, 5.20265, 0.0830948, 0.0481, 344.5086, 2446800.5, -9.4, 98.44, jup404, 0.0, 0.0, 0.0);

      internal static orbit saturn = new orbit("Saturn", 2446800.5, 2.4858, 113.555, 337.969, 9.5405, 0.033451, 0.052786, 159.6327, 2446800.5, -8.88, 82.73, sat404, 0.0, 0.0, 0.0);

      internal static orbit uranus = new orbit("Uranus", 2446800.5, 0.7738, 73.994, 98.746, 19.2233, 0.0116943, 0.045682, 84.8516, 2446800.5, -7.19, 35.02, ura404, 0.0, 0.0, 0.0);

      internal static orbit neptune = new orbit("Neptune", 2446800.5, 1.7697,
         131.677, 250.623, 30.1631,
         0.00594978,
         0.009019,
         254.2568, 2446800.5, -6.87,
         33.5, nep404, 0.0, 0.0, 0.0);

      internal static orbit pluto = new orbit("Pluto", 2446640.5, 17.1346, 110.204, 114.21, 39.4633, 0.0039757, 0.248662, 355.0554, 2446640.5, -1.0, 2.07, plu404, 0.0, 0.0, 0.0);
      #endregion

      #region GPlan
      private static double[,] ss = new double[NARGS + 1, 32];
      private static double[,] cc = new double[NARGS + 1, 32];
      private static double[] Args = new double[NARGS + 1];
      private static double LP_equinox;
      private static double NF_arcsec;
      private static double Ea_arcsec;
      private static double pA_precession;

      ///*   Routines to chew through tables of perturbations.  */
      internal static double mods3600(double x)
      {
         return ((x) - 1296000.0 * Math.Floor((x) / 1296000.0));
      }

      ///* From Simon et al (1994)  */
      ///* Arc sec per 10000 Julian years.  */
      internal static double[] freqs = new double[] { 53810162868.8982, 21066413643.3548, 12959774228.3429, 6890507749.3988, 1092566037.7991, 439960985.5372, 154248119.3933, 78655032.0744, 52272245.1795 };

      ///* Arc sec.  */
      internal static double[] phases = new double[] { 252.25090552 * 3600.0, 181.97980085 * 3600.0, 100.46645683 * 3600.0, 355.43299958 * 3600.0, 34.35151874 * 3600.0, 50.0774443 * 3600.0, 314.05500511 * 3600.0, 304.34866548 * 3600.0, 860492.1546 };

      internal static int gplan(double JD, plantbl plan, ref double[] pobj)
      {
         double su = 0;
         double cu = 0;
         double sv = 0;
         double cv = 0;
         double TI = 0;
         double t = 0;
         double sl = 0;
         double sb = 0;
         double sr = 0;
         int i = 0;
         int j = 0;
         int k = 0;
         int m = 0;
         int n = 0;
         int k1;
         int ip = 0;
         int np = 0;
         int nt = 0;
         int p = 0;
         int pl = 0;
         int pb = 0;
         int pr = 0;

         TI = System.Convert.ToDouble((JD - J2000) / plan.timescale);
         n = plan.maxargs;
         ///* Calculate sin( i*MM ), etc. for needed multiple angles.  */
         for (i = 0; i <= n - 1; i++) {
            j = plan.max_harmonic[i];
            if (j > 0) {
               sr = System.Convert.ToDouble((mods3600(freqs[i] * TI) + phases[i]) * STR);
               sscc(i, sr, j);
            }
         }

         ///* Point to start of table of arguments. */

         p = 0; //p = plan.arg_tbl
                ///* Point to tabulated cosine and sine amplitudes.  */
         pl = 0; //pl = plan.lon_tbl
         pb = 0; //pb = plan.lat_tbl
         pr = 0; //pr = plan.rad_tbl

         sl = 0.0;
         sb = 0.0;
         sr = 0.0;

         do {
            ///* argument of sine and cosine */
            ///* Number of periodic arguments. */
            np = plan.arg_tbl[p];
            p++;
            if (np < 0) {
               break;
            }
            if (np == 0) ///* It is a polynomial term.  */
        {
               nt = plan.arg_tbl[p];
               p++;
               cu = plan.lon_tbl[pl];
               pl++; ///* Longitude polynomial. */
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.lon_tbl[pl];
                  pl++;
               }
               sl += mods3600(cu);

               cu = plan.lat_tbl[pb];
               pb++; ///* Latitude polynomial. */
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.lat_tbl[pb];
                  pb++;
               }
               sb += cu;

               cu = plan.rad_tbl[pr];
               pr++; ///* Radius polynomial. */
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.rad_tbl[pr];
                  pr++;
               }
               sr += cu;
            }
            else {
               k1 = 0;
               cv = 0.0;
               sv = 0.0;
               for (ip = 0; ip <= np - 1; ip++) {
                  j = plan.arg_tbl[p];
                  p++; ///* What harmonic.  */
                  m = plan.arg_tbl[p] - 1;
                  p++; ///* Which planet.  */
                  if (j != 0) {
                     k = j;
                     if (j < 0) {
                        k = System.Convert.ToInt32(-k);
                     }
                     k--;
                     su = ss[m, k]; ///* sin(k*angle) */
                     if (j < 0) {
                        su = System.Convert.ToDouble(-su);
                     }
                     cu = cc[m, k];
                     if (k1 == 0) {
                        ///* set first angle */
                        sv = su;
                        cv = cu;
                        k1 = 1;
                     }
                     else {
                        ///* combine angles */
                        t = su * cv + cu * sv;
                        cv = cu * cv - su * sv;
                        sv = t;
                     }
                  }
               }
               ///* Highest power of T.  */
               nt = plan.arg_tbl[p];
               p++;
               cu = plan.lon_tbl[pl];
               pl++; ///* Longitude. */
               su = plan.lon_tbl[pl];
               pl++;
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.lon_tbl[pl];
                  pl++;
                  su = su * TI + plan.lon_tbl[pl];
                  pl++;
               }
               sl += cu * cv + su * sv;

               cu = plan.lat_tbl[pb];
               pb++; ///* Latitiude. */
               su = plan.lat_tbl[pb];
               pb++;
               for (ip = 1; ip <= nt; ip++) {
                  cu = cu * TI + plan.lat_tbl[pb];
                  pb++;
                  su = su * TI + plan.lat_tbl[pb];
                  pb++;
               }
               sb += cu * cv + su * sv;

               cu = plan.rad_tbl[pr];
               pr++; ///* Radius. */
               su = plan.rad_tbl[pr];
               pr++;
               for (ip = 1; ip <= nt; ip++) {
                  cu = cu * TI + plan.rad_tbl[pr];
                  pr++;
                  su = su * TI + plan.rad_tbl[pr];
                  pr++;
               }
               sr += cu * cv + su * sv;
            }
         } while (true);

         pobj[0] = STR * sl;
         pobj[1] = STR * sb;
         pobj[2] = STR * plan.distance * sr + plan.distance;

         return 0;
      }


      ///* Prepare lookup table of sin and cos ( i*Lj )
      // * for required multiple angles
      // */
      internal static int sscc(int k, double arg, int n)
      {
         double cu = 0;
         double su = 0;
         double cv = 0;
         double sv = 0;
         double s = 0;
         int i = 0;

         su = System.Convert.ToDouble(Math.Sin(arg));
         cu = System.Convert.ToDouble(Math.Cos(arg));
         ss[k, 0] = su; ///* sin(L) */
         cc[k, 0] = cu; ///* cos(L) */
         sv = System.Convert.ToDouble(2.0 * su * cu);
         cv = cu * cu - su * su;
         ss[k, 1] = sv; ///* sin(2L) */
         cc[k, 1] = cv;
         for (i = 2; i <= n - 1; i++) {
            s = su * cv + cu * sv;
            cv = cu * cv - su * sv;
            sv = s;
            ss[k, i] = sv; ///* sin( i+1 L ) */
            cc[k, i] = cv;
         }
         return 0;
      }
      ///* Compute mean elements at Julian date J.  */

      static public void mean_elements(double J)
      {
         double x = 0;
         double T = 0;
         double T2 = 0;

         ///* Time variables.  T is in Julian centuries.  */
         T = System.Convert.ToDouble((J - 2451545.0) / 36525.0);
         T2 = T * T;

         ///* Mean longitudes of planets (Simon et al, 1994) .047" subtracted from constant term for offset to DE403 origin. */

         ///* Mercury */
         x = mods3600(538101628.68898189 * T + 908103.213);
         x += System.Convert.ToDouble((0.00000639 * T - 0.0192789) * T2);
         Args[0] = STR * x;

         ///* Venus */
         x = mods3600(210664136.43354821 * T + 655127.236);
         x += System.Convert.ToDouble((-0.00000627 * T + 0.0059381) * T2);
         Args[1] = STR * x;

         ///* Earth  */
         x = mods3600(129597742.283429 * T + 361679.198);
         x += System.Convert.ToDouble((-0.00000523 * T - 0.0204411) * T2);
         Ea_arcsec = x;
         Args[2] = STR * x;

         ///* Mars */
         x = mods3600(68905077.493988 * T + 1279558.751);
         x += System.Convert.ToDouble((-0.00001043 * T + 0.0094264) * T2);
         Args[3] = STR * x;

         ///* Jupiter */
         x = mods3600(10925660.377991 * T + 123665.42);
         x += System.Convert.ToDouble(((((-0.00000000034 * T + 0.0000000591) * T + 0.000004667) * T + 0.00005706) * T - 0.3060378) * T2);
         Args[4] = STR * x;

         ///* Saturn */
         x = mods3600(4399609.855372 * T + 180278.752);
         x += System.Convert.ToDouble(((((0.00000000083 * T - 0.0000001452) * T - 0.000011484) * T - 0.00016618) * T + 0.7561614) * T2);
         Args[5] = STR * x;

         ///* Uranus */
         x = System.Convert.ToDouble(mods3600(1542481.193933 * T + 1130597.971) + (0.00002156 * T - 0.0175083) * T2);
         Args[6] = STR * x;

         ///* Neptune */
         x = System.Convert.ToDouble(mods3600(786550.320744 * T + 1095655.149) + (-0.00000895 * T + 0.0021103) * T2);
         Args[7] = STR * x;

         ///* Copied from cmoon.c, DE404 version.  */
         ///* Mean elongation of moon = D */
         x = mods3600(1602961600.9939659 * T + 1072261.2202445078);
         x += System.Convert.ToDouble((((((-0.0000000000003207663637426 * T + 0.00000000002555243317839) * T + 0.000000002560078201452) * T - 0.00003702060118571) * T + 0.0069492746836058421) * T - 6.7352202374457519) * T2); ///* D, t^2 */
         Args[9] = STR * x;

         ///* Mean distance of moon from its ascending node = F */
         x = mods3600(1739527262.8437717 * T + 335779.5141288474);
         x += System.Convert.ToDouble((((((0.0000000000004474984866301 * T + 0.00000000004189032191814) * T - 0.000000002790392351314) * T - 0.000002165750777942) * T - 0.00075311878482337989) * T - 13.117809789650071) * T2); ///* F, t^2 */
         NF_arcsec = x;
         Args[10] = STR * x;

         ///* Mean anomaly of sun = l' (J. Laskar) */
         x = mods3600(129596581.0230432 * T + 1287102.7407441526);
         x += System.Convert.ToDouble(((((((((1.62E-20 * T - 1.039E-17) * T - 0.00000000000000383508) * T + 0.0000000000004237343) * T + 0.000000000088555011) * T - 0.0000000477258489) * T - 0.000011297037031) * T + 0.0000874737173673247) * T - 0.55281306421783094) * T2);
         Args[11] = STR * x;

         ///* Mean anomaly of moon = l */
         x = mods3600(1717915922.8846793 * T + 485868.17465825332);
         x += System.Convert.ToDouble((((((-0.000000000001755312760154) * T + 0.00000000003452144225877 * T - 0.00000002506365935364) * T - 0.0002536291235258) * T + 0.052099641302735818) * T + 31.501359071894147) * T2); // /* l, t^2 */
         Args[12] = STR * x;

         ///* Mean longitude of moon, re mean ecliptic and equinox of date = L  */
         x = mods3600(1732564372.0442266 * T + 785939.8092105242);
         x += System.Convert.ToDouble((((((0.00000000000007200592540556 * T + 0.0000000002235210987108) * T - 0.00000001024222633731) * T - 0.00006073960534117) * T + 0.006901724852838049) * T - 5.65504600274714) * T2); // /* L, t^2 */
         LP_equinox = x;
         Args[13] = STR * x;

         ///* Precession of the equinox  */
         x = System.Convert.ToDouble((((((((((-8.66E-20 * T - 4.759E-17) * T + 0.000000000000002424) * T + 0.0000000000013095) * T + 0.00000000017451) * T - 0.000000018055) * T - 0.0000235316) * T + 0.000076) * T + 1.105414) * T + 5028.791959) * T);
         ///* Moon's longitude re fixed J2000 equinox.  */
         ///*
         //Args(13) -= x;
         //*/
         pA_precession = STR * x;

         ///* Free librations.  */
         ///* longitudinal libration 2.891725 years */
         x = mods3600(44817540.9 * T + 806045.7);
         Args[14] = STR * x;
         ///* libration P, 24.2 years */
         x = mods3600(5364867.87 * T - 391702.8);
         Args[15] = STR * x;

         //Args(16) = 0.0

         ///* libration W, 74.7 years. */
         x = mods3600(1735730.0 * T);
         Args[17] = STR * x;
      }


      ///* Generic program to accumulate sum of trigonometric series
      //   in three variables (e.g., longitude, latitude, radius)
      //   of the same list of arguments.  */

      internal static int g3plan(double JD, plantbl plan, ref double[] pobj, int objnum)
      {
         int i = 0;
         int j = 0;
         int k = 0;
         int m = 0;
         int n = 0;
         int k1;
         int ip = 0;
         int np = 0;
         int nt = 0;
         int p = 0;
         int pl = 0;
         int pb = 0;
         int pr = 0;
         double su = 0;
         double cu = 0;
         double sv = 0;
         double cv = 0;
         double TI = 0;
         double t = 0;
         double sl = 0;
         double sb = 0;
         double sr = 0;

         mean_elements(JD);
         //#If 0 Then
         //  /* For librations, moon's longitude is sidereal.  */
         //            If (flag) Then
         //    Args(13) -= pA_precession;
         //#End If

         TI = System.Convert.ToDouble((JD - J2000) / plan.timescale);
         n = plan.maxargs;
         ///* Calculate sin( i*MM ), etc. for needed multiple angles.  */
         for (i = 0; i <= n - 1; i++) {
            j = plan.max_harmonic[i];
            if (j > 0) {
               sscc(i, Args[i], j);
            }
         }

         ///* Point to start of table of arguments. */
         p = 0; //plan.arg_tbl
                ///* Point to tabulated cosine and sine amplitudes.  */
         pl = 0; //plan.lon_tbl
         pb = 0; //plan.lat_tbl
         pr = 0; //plan.rad_tbl
         sl = 0.0;
         sb = 0.0;
         sr = 0.0;

         do {
            ///* argument of sine and cosine */
            ///* Number of periodic arguments. */
            np = plan.arg_tbl[p];
            p++;
            if (np < 0) {
               break;
            }
            if (np == 0) ///* It is a polynomial term.  */
        {
               nt = plan.arg_tbl[p];
               p++;
               cu = plan.lon_tbl[pl];
               pl++; ///* "Longitude" polynomial (phi). */
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.lon_tbl[pl];
                  pl++;
               }
               ///*	  sl +=  mods3600 (cu); */
               sl += cu;

               cu = plan.lat_tbl[pb];
               pb++; ///* "Latitude" polynomial (theta). */
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.lat_tbl[pb];
                  pb++;
               }
               sb += cu;

               cu = plan.rad_tbl[pr];
               pr++; ///* Radius polynomial (psi). */
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.rad_tbl[pr];
                  pr++;
               }
               sr += cu;
            }
            else {
               k1 = 0;
               cv = 0.0;
               sv = 0.0;
               for (ip = 0; ip <= np - 1; ip++) {
                  j = plan.arg_tbl[p];
                  p++; ///* What harmonic.  */
                  m = plan.arg_tbl[p] - 1;
                  p++; ///* Which planet.  */
                  if (j != 0) {
                     ///*	      k = abs (j); */
                     if (j < 0) {
                        k = System.Convert.ToInt32(-j);
                     }
                     else {
                        k = j;
                     }
                     k--;
                     su = ss[m, k]; ///* sin(k*angle) */
                     if (j < 0) {
                        su = System.Convert.ToDouble(-su);
                     }
                     cu = cc[m, k];
                     if (k1 == 0) ///* set first angle */
                 {
                        sv = su;
                        cv = cu;
                        k1 = 1;
                     }
                     else {
                        ///* combine angles */
                        t = su * cv + cu * sv;
                        cv = cu * cv - su * sv;
                        sv = t;
                     }
                  }
               }
               ///* Highest power of T.  */
               nt = plan.arg_tbl[p];
               p++;

               ///* Longitude. */
               cu = plan.lon_tbl[pl];
               pl++;
               su = plan.lon_tbl[pl];
               pl++;
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.lon_tbl[pl];
                  pl++;
                  su = su * TI + plan.lon_tbl[pl];
                  pl++;
               }
               sl += cu * cv + su * sv;

               ///* Latitiude. */
               cu = plan.lat_tbl[pb];
               pb++;
               su = plan.lat_tbl[pb];
               pb++;
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.lat_tbl[pb];
                  pb++;
                  su = su * TI + plan.lat_tbl[pb];
                  pb++;
               }
               sb += cu * cv + su * sv;

               ///* Radius. */
               cu = plan.rad_tbl[pr];
               pr++;
               su = plan.rad_tbl[pr];
               pr++;
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.rad_tbl[pr];
                  pr++;
                  su = su * TI + plan.rad_tbl[pr];
                  pr++;
               }
               sr += cu * cv + su * sv;
            }
         } while (true);
         t = plan.trunclvl;
         pobj[0] = System.Convert.ToDouble(Args[objnum - 1] + STR * t * sl);
         pobj[1] = STR * t * sb;
         pobj[2] = plan.distance * (1.0 + STR * t * sr);
         return 0;
      }

      ///* Generic program to accumulate sum of trigonometric series
      //   in two variables (e.g., longitude, radius)
      //   of the same list of arguments.  */
      internal static int g2plan(double JD, plantbl plan, ref double[] pobj)
      {
         int i = 0;
         int j = 0;
         int k = 0;
         int m = 0;
         int n = 0;
         int k1;
         int ip = 0;
         int np = 0;
         int nt = 0;
         int p = 0;
         int pl = 0;
         int pr = 0;
         double su = 0;
         double cu = 0;
         double sv = 0;
         double cv = 0;
         double TI = 0;
         double t = 0;
         double sl = 0;
         double sr = 0;

         mean_elements(JD);
         //#If 0 Then
         ///* For librations, moon's longitude is sidereal.  */
         //If (flag) Then
         //Args(13) -= pA_precession;
         //#End If
         TI = System.Convert.ToDouble((JD - J2000) / plan.timescale);
         n = plan.maxargs;
         ///* Calculate sin( i*MM ), etc. for needed multiple angles.  */
         for (i = 0; i <= n - 1; i++) {
            j = plan.max_harmonic[i];
            if (j > 0) {
               sscc(i, Args[i], j);
            }
         }

         ///* Point to start of table of arguments. */
         p = 0; //plan.arg_tbl
                ///* Point to tabulated cosine and sine amplitudes.  */
         pl = 0; //(long *) plan.lon_tbl;
         pr = 0; //(long *) plan.rad_tbl;
         sl = 0.0;
         sr = 0.0;

         do {
            ///* argument of sine and cosine */
            ///* Number of periodic arguments. */
            np = plan.arg_tbl[p];
            p++; //*p++;
            if (np < 0) {
               break;
            }

            if (np == 0) ///* It is a polynomial term.  */
        {
               nt = plan.arg_tbl[p];
               p++;
               cu = plan.lon_tbl[pl];
               pl++; //*pl++; '/* Longitude polynomial. */
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.lon_tbl[pl];
                  pl++; //*pl++;
               }
               ///*	  sl +=  mods3600 (cu); */
               sl += cu;
               ///* Radius polynomial. */
               cu = plan.rad_tbl[pr];
               pr++; //*pr++;
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.rad_tbl[pr];
                  pr++;
               }
               sr += cu;
            }
            else {
               k1 = 0;
               cv = 0.0;
               sv = 0.0;
               for (ip = 0; ip <= np - 1; ip++) {
                  j = plan.arg_tbl[p];
                  p++; ///* What harmonic.  */
                  m = plan.arg_tbl[p] - 1;
                  p++; ///* Which planet.  */
                  if (j != 0) {
                     ///*	      k = abs (j); */
                     if (j < 0) {
                        k = System.Convert.ToInt32(-j);
                     }
                     else {
                        k = j;
                     }
                     k--;
                     su = ss[m, k]; ///* sin(k*angle) */
                     if (j < 0) {
                        su = System.Convert.ToDouble(-su);
                     }
                     cu = cc[m, k];
                     if (k1 == 0) {
                        ///* set first angle */
                        sv = su;
                        cv = cu;
                        k1 = 1;
                     }
                     else {
                        ///* combine angles */
                        t = su * cv + cu * sv;
                        cv = cu * cv - su * sv;
                        sv = t;
                     }
                  }
               }
               ///* Highest power of T.  */
               nt = plan.arg_tbl[p];
               p++; //*p++;
                    ///* Longitude. */
               cu = plan.lon_tbl[pl];
               pl++;
               su = plan.lon_tbl[pl];
               pl++;
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.lon_tbl[pl];
                  pl++;
                  su = su * TI + plan.lon_tbl[pl];
                  pl++;
               }
               sl += cu * cv + su * sv;
               ///* Radius. */
               cu = plan.rad_tbl[pr];
               pr++;
               su = plan.rad_tbl[pr];
               pr++;
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.rad_tbl[pr];
                  pr++;
                  su = su * TI + plan.rad_tbl[pr];
                  pr++;
               }
               sr += cu * cv + su * sv;
            }
         } while (true);
         t = plan.trunclvl;
         pobj[0] = t * sl;
         pobj[2] = t * sr;
         return 0;
      }


      ///* Generic program to accumulate sum of trigonometric series
      //   in one variable.  */

      internal static double g1plan(double JD, plantbl plan)
      {
         int i = 0;
         int j = 0;
         int k = 0;
         int m = 0;
         int k1;
         int ip = 0;
         int np = 0;
         int nt = 0;
         int p = 0;
         int pl = 0;
         double su = 0;
         double cu = 0;
         double sv = 0;
         double cv = 0;
         double TI = 0;
         double t = 0;
         double sl = 0;

         TI = System.Convert.ToDouble((JD - J2000) / plan.timescale);
         mean_elements(JD);
         ///* Calculate sin( i*MM ), etc. for needed multiple angles.  */
         for (i = 0; i <= NARGS - 1; i++) {
            j = plan.max_harmonic[i];
            if (j > 0) {
               sscc(i, Args[i], j);
            }
         }

         ///* Point to start of table of arguments. */
         p = 0; //plan.arg_tbl;
                ///* Point to tabulated cosine and sine amplitudes.  */
         pl = 0; //(long *) plan.lon_tbl;
         sl = 0.0;

         do ///* argument of sine and cosine */
      {
            ///* Number of periodic arguments. */
            np = plan.arg_tbl[p];
            p++; //*p++;
            if (np < 0) {
               break;
            }
            if (np == 0) {
               ///* It is a polynomial term.  */
               nt = plan.arg_tbl[p];
               p++; //*p++;
               cu = plan.lon_tbl[pl];
               pl++; //*pl++;
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.lon_tbl[pl];
                  pl++; //*pl++;
               }
               ///*	  sl +=  mods3600 (cu); */
               sl += cu;
            }
            else {
               k1 = 0;
               cv = 0.0;
               sv = 0.0;
               for (ip = 0; ip <= np - 1; ip++) {
                  j = plan.arg_tbl[p];
                  p++; ///* What harmonic.  */
                  m = plan.arg_tbl[p] - 1;
                  p++; ///* Which planet.  */
                  if (j != 0) {
                     ///*	      k = abs (j); */
                     if (j < 0) {
                        k = System.Convert.ToInt32(-j);
                     }
                     else {
                        k = j;
                     }
                     k--;
                     su = ss[m, k]; ///* sin(k*angle) */
                     if (j < 0) {
                        su = System.Convert.ToDouble(-su);
                     }
                     cu = cc[m, k];
                     if (k1 == 0) {
                        ///* set first angle */
                        sv = su;
                        cv = cu;
                        k1 = 1;
                     }
                     else {
                        ///* combine angles */
                        t = su * cv + cu * sv;
                        cv = cu * cv - su * sv;
                        sv = t;
                     }
                  }
               }
               ///* Highest power of T.  */
               nt = plan.arg_tbl[p];
               p++;
               ///* Math.Cosine and sine coefficients.  */
               cu = plan.lon_tbl[pl];
               pl++;
               su = plan.lon_tbl[pl];
               pl++;
               for (ip = 0; ip <= nt - 1; ip++) {
                  cu = cu * TI + plan.lon_tbl[pl];
                  pl++;
                  su = su * TI + plan.lon_tbl[pl];
                  pl++;
               }
               sl += cu * cv + su * sv;
            }
         } while (true);
         return (plan.trunclvl * sl);
      }

      internal static int gmoon(double J, ref double[] rect, ref double[] pol)
      {
         double x = 0;
         double cosB = 0;
         double sinB = 0;
         double cosL = 0;
         double sinL = 0;
         double eps = 0;
         double coseps = 0;
         double sineps = 0;
         //Dim TL As New TraceLogger("", "GMoon")
         //TL.Enabled = True
         //TL.LogMessage("GMoon", "Before G2Plan")
         g2plan(J, moonlr, ref pol);
         //TL.LogMessage("GMoon", "After G2Plan")
         x = pol[0];
         x += LP_equinox;
         if (x < -645000.0) {
            x += 1296000.0;
         }
         if (x > 645000.0) {
            x -= 1296000.0;
         }
         pol[0] = STR * x;
         //TL.LogMessage("GMoon", "Before G1Plan")
         x = g1plan(J, moonlat);
         //TL.LogMessage("GMoon", "After G1Plan")
         pol[1] = STR * x;
         x = System.Convert.ToDouble((1.0 + STR * pol[2]) * moonlr.distance);
         pol[2] = x;
         ///* Convert ecliptic polar to equatorial rectangular coordinates.  */
         //TL.LogMessage("GMoon", "Before Epsilin")
         epsiln(J, ref eps, ref coseps, ref sineps);
         //TL.LogMessage("GMoon", "After Epsilin")
         cosB = System.Convert.ToDouble(Math.Cos(pol[1]));
         sinB = System.Convert.ToDouble(Math.Sin(pol[1]));
         cosL = System.Convert.ToDouble(Math.Cos(pol[0]));
         sinL = System.Convert.ToDouble(Math.Sin(pol[0]));
         rect[0] = cosB * cosL * x;
         rect[1] = System.Convert.ToDouble((coseps * cosB * sinL - sineps * sinB) * x);
         rect[2] = System.Convert.ToDouble((sineps * cosB * sinL + coseps * sinB) * x);
         //TL.Enabled = False : TL.Dispose() : TL = Nothing

         return 0;
      }

      #endregion
   }
}