using System;

namespace ASCOM.Astrometry
{
   sealed class EphemerisCode
   {

      //Function patterned after get_earth() in the original NOVAS-C V2 package.
      //This function returns (via ref-params) the barycentric TDT and both
      //heliocentric and barycentric position and velocity of the Earth, at
      //the given TJD. You can pass an IDispatch pointer for an ephemeris
      //component, and it will be used. If that is NULL the internal solsys3()
      //function is used (see solarsystem().
      //
      //For more info, see the original NOVAS-C sources.
      // VBConversions Note: Former VB static variables moved to class level because they aren't supported in C#.
      static double get_earth_nov_tjd_last = 0.0;

      internal static void get_earth_nov(ref IEphemeris pEphDisp, double tjd, ref double tdb, ref double[] peb, ref double[] veb, ref double[] pes, ref double[] ves)
      {
         short i = 0;
         short rc = 0;
         double dummy = 0;
         double secdiff = 0;
         // static double tjd_last = 0.0; VBConversions Note: Static variable moved to class level and renamed get_earth_nov_tjd_last. Local static variables are not supported in C#.
         double ltdb = 0;
         double[] lpeb = new double[4];
         double[] lveb = new double[4];
         double[] lpes = new double[4];
         double[] lves = new double[4];
         //Dim TL As New TraceLogger("", "get_earth_nov")
         //TL.Enabled = True
         //TL.LogMessage("get_earth_nov", "Start")
         ////
         //// Compute the TDB Julian date corresponding to 'tjd'.
         ////

         //If (Abs(tjd - tjd_last) > 0.000001) Then 'Optimize repeated calls
         Tdb2Tdt(tjd, dummy, secdiff);
         //TL.LogMessage("get_earth_nov", "after tbd2tdt")
         ltdb = tjd + secdiff / 86400.0;

         ////
         //// Get position and velocity of the Earth wrt barycenter of
         //// solar system and wrt center of the sun. These calls reflect
         //// exceptions thrown by the attached ephemeris generator, so
         //// we just return the hr ... the ErrorInfo is already set!
         ////
         try {
            //TL.LogMessage("get_earth_nov", "before solsysnov barycentric")
            rc = solarsystem_nov(ref pEphDisp, tjd, ltdb, Body.Earth, Origin.Barycentric, ref lpeb, ref lveb);
            //TL.LogMessage("get_earth_nov", "after solsysnov barycentric")
            if (rc != 0) {
               throw (new Exceptions.NOVASFunctionException("EphemerisCode:get_earth_nov Earth eph exception", "solarsystem_nov", rc));
            }
         }
         catch (Exception) {
            get_earth_nov_tjd_last = 0.0;
            throw;
         }

         try {
            //TL.LogMessage("get_earth_nov", "before solsysnov heliocentric")
            rc = solarsystem_nov(ref pEphDisp, tjd, ltdb, Body.Earth, Origin.Heliocentric, ref lpes, ref lves);
            //TL.LogMessage("get_earth_nov", "after solsysnov heliocentric")
            if (rc != 0) {
               throw (new Exceptions.NOVASFunctionException("EphemerisCode:get_earth_nov Earth eph exception", "solarsystem_nov", rc));
            }
         }
         catch (Exception) {
            get_earth_nov_tjd_last = 0.0;
            throw;
         }

         get_earth_nov_tjd_last = tjd;
         //End If
         tdb = ltdb;
         for (i = 0; i <= 2; i++) {
            peb[i] = lpeb[i];
            veb[i] = lveb[i];
            pes[i] = lpes[i];
            ves[i] = lves[i];
         }

         //TL.Enabled = False
         //TL.Dispose()
         //TL = Nothing

      }

      ////
      //// Ephemeris() - Wrapper for external ephemeris generator
      ////
      //// The ephemeris generator must support a single method:
      ////
      ////     result(6) = GetPositionAndVelocity(tjd, Type, Number, Name)
      ////
      ////	tjd		Terrestrial Julian Date
      ////	Type	Type of body: 0 = major planet, Sun, or Moon
      ////						  1 = minor planet
      ////	Number: For Type = 0: Mercury = 1, ..., Pluto = 9
      ////			For Type = 1: minor planet number or 0 for unnumbered MP
      ////  Name:   For Type = 0: n/a
      ////			For Type = 1: n/a for numbered MPs. For unnumbered MPs, this
      ////						  is the MPC PACKED designation.
      ////  result	A SAFEARRAY of VARIANT, each element VT_R8 (double). Elements
      ////			0-2 are the position vector of the body, elements 3.5 are the
      ////			velocity vector of the body.
      ////
      internal static void ephemeris_nov(IEphemeris ephDisp, double tjd, BodyType btype, int num, string name, Origin origin, ref double[] pos, ref double[] vel)
      {
         int i = 0;
         double[] posvel = new double[7];
         double[] p = new double[3];
         double[] v = new double[3];
         //Dim bdy As bodystruct
         //Dim org As NOVAS2Net.Origin
         //Dim rc As Short
         //Dim TL As New TraceLogger("", "EphNov")
         //TL.Enabled = True
         // TL.LogMessage("EphNov", "Start")
         ////
         //// Check inputs
         ////
         if (ephDisp == null) {
            throw (new Exceptions.ValueNotSetException("Ephemeris_nov Ephemeris object not set"));
         }
         else {
            if ((origin != origin.Barycentric) && (origin != origin.Heliocentric)) {
               throw (new Utilities.Exceptions.InvalidValueException("Ephemeris_nov Origin is neither barycentric or heliocentric"));
            }

            ////
            //// Call the ephemeris for the heliocentric J2000.0 equatorial coordinates
            BodyType kbtype = default(BodyType);
            //TL.LogMessage("EphNov", "Before Case Btype")
            if (btype == BodyType.Comet) {
               kbtype = BodyType.Comet;
            }
            else if (btype == BodyType.MajorPlanet) {
               kbtype = BodyType.MajorPlanet;
            }
            else if (btype == BodyType.MinorPlanet) {
               kbtype = BodyType.MinorPlanet;
            }

            Body knum = default(Body);
            switch (num) {
               case 1:
                  knum = Body.Mercury;
                  break;
               case 2:
                  knum = Body.Venus;
                  break;
               case 3:
                  knum = Body.Earth;
                  break;
               case 4:
                  knum = Body.Mars;
                  break;
               case 5:
                  knum = Body.Jupiter;
                  break;
               case 6:
                  knum = Body.Saturn;
                  break;
               case 7:
                  knum = Body.Uranus;
                  break;
               case 8:
                  knum = Body.Neptune;
                  break;
               case 9:
                  knum = Body.Pluto;
                  break;
            }
            ephDisp.BodyType = kbtype;
            ephDisp.Number = knum;
            if (name != "") {
               ephDisp.Name = name;
            }
            //TL.LogMessage("EphNov", "Before ephDisp GetPosAndVel")
            posvel = ephDisp.GetPositionAndVelocity(tjd);
            //TL.LogMessage("EphNov", "After ephDisp GetPosAndVel")
         }

         if (origin == origin.Barycentric) {

            double[] sun_pos = new double[4];
            double[] sun_vel = new double[4];

            //// CHICKEN AND EGG ALERT!!! WE CANNOT CALL OURSELVES FOR
            //// BARYCENTER CALCULATION -- AS AN APPROXIMATION, WE USE
            //// OUR INTERNAL SOLSYS3() FUNCTION TO GET THE BARYCENTRIC
            //// SUN. THIS SHOULD BE "GOOD ENOUGH". IF WE EVER GET
            //// AN EPHEMERIS GEN THAT HANDLES BARYCENTRIC, WE CAN
            //// CAN THIS...
            //TL.LogMessage("EphNov", "Before solsys3")
            solsys3_nov(tjd, Body.Sun, origin.Barycentric, ref sun_pos, ref sun_vel);
            //TL.LogMessage("EphNov", "After solsys3")
            for (i = 0; i <= 2; i++) {
               posvel[i] += sun_pos[i];
               posvel[i + 3] += (int)(sun_vel[i]);
            }
         }

         for (i = 0; i <= 2; i++) {
            pos[i] = posvel[i];
            vel[i] = posvel[i + 3];
         }
         //TL.Enabled = False
         //TL.Dispose()
      }

      //// ===============
      //// LOCAL FUNCTIONS
      //// ===============


      ////
      //// This is the function used to get the position and velocity vectors
      //// for the major solar system bodies and the moon. It is patterned after
      //// the solarsystem() function in the original NOVAS-C package. You can
      //// pass an IDispatch pointer for an ephemeris component, and it will be
      //// used. If that is NULL the internal solsys3() function is used.
      ////
      //// This function must set error info... it is designed to work with
      //// reflected exceptions from the attached ephemeris
      ////
      internal static short solarsystem_nov(ref IEphemeris ephDisp, double tjd, double tdb, Body planet, Origin origin, ref double[] pos, ref double[] vel)
      {
         //Dim pl As NOVAS2.Body, org As NOVAS2.Origin
         //Dim TL As New TraceLogger("", "solarsystem_nov")
         short rc = 0;
         //TL.Enabled = True
         //TL.LogMessage("solarsystem_nov", "Start")
         ////
         //// solsys3 takes tdb, ephemeris takes tjd
         ////
         //Select Case origin
         //    Case OriginType.nvBarycentric
         //org = NOVAS2.Origin.SolarSystemBarycentre
         //    Case OriginType.nvHeliocentric
         //org = NOVAS2.Origin.CentreOfMassOfSun
         //End Select
         //Select Case planet
         //    Case PlanetNumber.nvEarth
         //pl = NOVAS2.Body.Earth
         //    Case PlanetNumber.nvJupiter
         //pl = NOVAS2.Body.Jupiter
         //    Case PlanetNumber.nvMars
         //pl = NOVAS2.Body.Mars
         //    Case PlanetNumber.nvMercury
         //pl = NOVAS2.Body.Mercury
         //    Case PlanetNumber.nvMoon
         //pl = NOVAS2.Body.Moon
         //    Case PlanetNumber.nvNeptune
         //pl = NOVAS2.Body.Neptune
         //    Case PlanetNumber.nvPluto
         //pl = NOVAS2.Body.Pluto
         //    Case PlanetNumber.nvSaturn
         //pl = NOVAS2.Body.Saturn
         //    Case PlanetNumber.nvSun
         //pl = NOVAS2.Body.Sun
         //    Case PlanetNumber.nvUranus
         //pl = NOVAS2.Body.Uranus
         //    Case PlanetNumber.nvVenus
         //pl = NOVAS2.Body.Venus
         //End Select
         //TL.LogMessage("solarsystem_nov", "After planet")
         if (ephDisp == null) //No ephemeris attached
         {
            //rc = solsys3_nov(tdb, planet, origin, pos, vel)
            throw (new Exceptions.ValueNotSetException("EphemerisCode:SolarSystem_Nov No emphemeris object supplied"));
         }
         else {
            //CHECK TDB BELOW IS CORRECT!
            //TL.LogMessage("solarsystem_nov", "Before ephemeris_nov")
            ephemeris_nov(ephDisp, tdb, BodyType.MajorPlanet, planet, "", origin, ref pos, ref vel);
            //TL.LogMessage("solarsystem_nov", "After ephemeris_nov")
         }
         //TL.Enabled = False
         //TL.Dispose()
         //TL = Nothing
         return rc;
      }

      ////
      //// solsys3() - Internal function that gives reasonable ephemerides for
      //// Sun or Earth, barycentric or heliocentric.
      ////
      // VBConversions Note: Former VB static variables moved to class level because they aren't supported in C#.
      static double solsys3_nov_tlast = 0.0;
      static double solsys3_nov_sine = 0;
      static double solsys3_nov_cose = 0;
      static double solsys3_nov_tmass = 0;
      static double[] solsys3_nov_pbary = new double[4];
      static double[] solsys3_nov_vbary = new double[4];

      private static short solsys3_nov(double tjd, Body body, Origin origin, ref double[] pos, ref double[] vel)
      {

         int i = 0;

         ///*
         //The arrays below contain data for the four largest planets.  Masses
         //are DE405 values; elements are from Explanatory Supplement, p. 316).
         //These data are used for barycenter computations only.
         //*/

         double[] pm = new double[] { 1047.349, 3497.898, 22903.0, 19412.2 };
         double[] pa = new double[] { 5.203363, 9.53707, 19.191264, 30.068963 };
         double[] pl = new double[] { 0.60047, 0.871693, 5.466933, 5.32116 };
         double[] pn = new double[] { 0.001450138, 0.0005841727, 0.0002047497, 0.0001043891 };

         ///*
         //obl' is the obliquity of ecliptic at epoch J2000.0 in degrees.
         //*/

         const double obl = 23.43929111;

         // static double tlast = 0.0; VBConversions Note: Static variable moved to class level and renamed solsys3_nov_tlast. Local static variables are not supported in C#.
         // static double sine = 0; VBConversions Note: Static variable moved to class level and renamed solsys3_nov_sine. Local static variables are not supported in C#.
         // static double cose = 0; VBConversions Note: Static variable moved to class level and renamed solsys3_nov_cose. Local static variables are not supported in C#.
         // static double tmass = 0; VBConversions Note: Static variable moved to class level and renamed solsys3_nov_tmass. Local static variables are not supported in C#.
         // static double[] pbary = new double[4]; VBConversions Note: Static variable moved to class level and renamed solsys3_nov_pbary. Local static variables are not supported in C#.
         // static double[] vbary = new double[4]; VBConversions Note: Static variable moved to class level and renamed solsys3_nov_vbary. Local static variables are not supported in C#.

         double oblr = 0;
         double qjd = 0;
         double ras = 0;
         double decs = 0;
         double diss = 0;
         double[] pos1 = new double[4];
         double[,] p = new double[4, 4];
         double dlon = 0;
         double sinl = 0;
         double cosl = 0;
         double x = 0;
         double y = 0;
         double z = 0;
         double xdot = 0;
         double ydot = 0;
         double zdot = 0;
         double f = 0;

         ////
         //// Check inputs
         ////
         if ((origin != origin.Barycentric) && (origin != origin.Heliocentric)) {
            throw (new Utilities.Exceptions.InvalidValueException("EphemerisCode.Solsys3 Invalid origin: " + System.Convert.ToString(origin)));
         }

         if ((tjd < 2340000.5) || (tjd > 2560000.5)) {
            throw (new Utilities.Exceptions.InvalidValueException("EphemerisCode.Solsys3 Invalid tjd: " + System.Convert.ToString(tjd)));
         }


         ///*
         //Initialize constants.
         //*/

         if (solsys3_nov_tlast == 0.0) {
            oblr = obl * TWOPI / 360.0;
            solsys3_nov_sine = System.Convert.ToDouble(Sin(oblr));
            solsys3_nov_cose = System.Convert.ToDouble(Cos(oblr));
            solsys3_nov_tmass = 1.0;
            for (i = 0; i <= 3; i++) {
               solsys3_nov_tmass += 1.0 / pm[i];
            }
            solsys3_nov_tlast = 1.0;
         }
         ///*
         //Form helicentric coordinates of the Sun or Earth, depending on
         //body'.
         //*/

         if ((body == 0) || (body == 1) || (body == 10)) {
            for (i = 0; i <= 2; i++) {
               pos[i] = 0.0;
               vel[i] = 0.0;
            }
         }
         else if ((body == 2) || (body == 3)) {
            for (i = 0; i <= 2; i++) {
               qjd = tjd + ((i) - 1.0) * 0.1;
               sun_eph_nov(qjd, ras, decs, diss);
               RADec2Vector(ras, decs, diss, pos1);
               Precession(qjd, pos1, T0, pos);
               p[i, 0] = -pos[0];
               p[i, 1] = -pos[1];
               p[i, 2] = -pos[2];
            }
            for (i = 0; i <= 2; i++) {
               pos[i] = p[1, i];
               vel[i] = System.Convert.ToDouble((p[2, i] - p[0, i]) / 0.2);
            }
         }
         else {
            throw (new Utilities.Exceptions.InvalidValueException("EphemerisCode.Solsys3 Invalid body: " + System.Convert.ToString(body)));
         }

         ///*
         //If 'origin' = 0, move origin to solar system barycenter.
         //
         //Solar system barycenter coordinates are computed from rough
         //approximations of the coordinates of the four largest planets.
         //*/

         if (origin == origin.Barycentric) {
            if (tjd != solsys3_nov_tlast) {
               for (i = 0; i <= 2; i++) {
                  solsys3_nov_pbary[i] = 0.0;
                  solsys3_nov_vbary[i] = 0.0;
               }

               ///*
               //The following loop cycles once for each of the four planets.
               //
               //sinl' and 'cosl' are the sine and cosine of the planet's mean
               //longitude.
               //*/

               for (i = 0; i <= 3; i++) {
                  dlon = System.Convert.ToDouble(pl[i] + pn[i] * (tjd - T0));
                  dlon = dlon % TWOPI;
                  sinl = System.Convert.ToDouble(Sin(dlon));
                  cosl = System.Convert.ToDouble(Cos(dlon));

                  x = System.Convert.ToDouble(pa[i] * cosl);
                  y = System.Convert.ToDouble(pa[i] * sinl * solsys3_nov_cose);
                  z = System.Convert.ToDouble(pa[i] * sinl * solsys3_nov_sine);
                  xdot = System.Convert.ToDouble(-pa[i] * pn[i] * sinl);
                  ydot = System.Convert.ToDouble(pa[i] * pn[i] * cosl * solsys3_nov_cose);
                  zdot = System.Convert.ToDouble(pa[i] * pn[i] * cosl * solsys3_nov_sine);

                  f = System.Convert.ToDouble(1.0 / (pm[i] * solsys3_nov_tmass));

                  solsys3_nov_pbary[0] += x * f;
                  solsys3_nov_pbary[1] += y * f;
                  solsys3_nov_pbary[2] += z * f;
                  solsys3_nov_vbary[0] += xdot * f;
                  solsys3_nov_vbary[1] += ydot * f;
                  solsys3_nov_vbary[2] += zdot * f;
               }

               solsys3_nov_tlast = tjd;
            }

            for (i = 0; i <= 2; i++) {
               pos[i] -= solsys3_nov_pbary[i];
               vel[i] -= solsys3_nov_vbary[i];
            }
         }
         return (short)0;
      }

      private struct sun_con
      {
         internal double l;
         internal double r;
         internal double alpha;
         internal double nu;
         internal sun_con(double pl, double pr, double palpha, double pnu)
         {
            l = pl;
            r = pr;
            alpha = palpha;
            nu = pnu;
         }
      }

      private static void sun_eph_nov(double jd, double ra, double dec, double dis)
      {
         int i = 0;

         double sum_lon = 0.0;
         double sum_r = 0.0;
         const double factor = 0.0000001;
         double u = 0;
         double arg = 0;
         double lon = 0;
         double lat;
         double t = 0;
         double t2 = 0;
         double emean = 0;
         double sin_lon = 0;

         sun_con[] con = new sun_con[] { new sun_con(403406.0, 0.0, 4.721964, 1.621043), new sun_con(195207.0, -97597.0, 5.937458, 62830.348067), new sun_con(119433.0, -59715.0, 1.115589, 62830.821524), new sun_con(112392.0, -56188.0, 5.781616, 62829.634302), new sun_con(3891.0, -1556.0, 5.5474, 125660.5691), new sun_con(2819.0, -1126.0, 1.512, 125660.9845), new sun_con(1721.0, -861.0, 4.1897, 62832.4766), new sun_con(0.0, 941.0, 1.163, 0.813), new sun_con(660.0, -264.0, 5.415, 125659.31), new sun_con(350.0, -163.0, 4.315, 57533.85), new sun_con(334.0, 0.0, 4.553, -33.931), new sun_con(314.0, 309.0, 5.198, 777137.715), new sun_con(268.0, -158.0, 5.989, 78604.191), new sun_con(242.0, 0.0, 2.911, 5.412), new sun_con(234.0, -54.0, 1.423, 39302.098), new sun_con(158.0, 0.0, 0.061, -34.861), new sun_con(132.0, -93.0, 2.317, 115067.698), new sun_con(129.0, -20.0, 3.193, 15774.337), new sun_con(114.0, 0.0, 2.828, 5296.67), new sun_con(99.0, -47.0, 0.52, 58849.27), new sun_con(93.0, 0.0, 4.65, 5296.11), new sun_con(86.0, 0.0, 4.35, -3980.7), new sun_con(78.0, -33.0, 2.75, 52237.69), new sun_con(72.0, -32.0, 4.5, 55076.47), new sun_con(68.0, 0.0, 3.23, 261.08), new sun_con(64.0, -10.0, 1.22, 15773.85), new sun_con(46.0, -16.0, 0.14, 188491.03), new sun_con(38.0, 0.0, 3.44, -7756.55), new sun_con(37.0, 0.0, 4.37, 264.89), new sun_con(32.0, -24.0, 1.14, 117906.27), new sun_con(29.0, -13.0, 2.84, 55075.75), new sun_con(28.0, 0.0, 5.96, -7961.39), new sun_con(27.0, -9.0, 5.09, 188489.81), new sun_con(27.0, 0.0, 1.72, 2132.19), new sun_con(25.0, -17.0, 2.56, 109771.03), new sun_con(24.0, -11.0, 1.92, 54868.56), new sun_con(21.0, 0.0, 0.09, 25443.93), new sun_con(21.0, 31.0, 5.98, -55731.43), new sun_con(20.0, -10.0, 4.03, 60697.74), new sun_con(18.0, 0.0, 4.27, 2132.79), new sun_con(17.0, -12.0, 0.79, 109771.63), new sun_con(14.0, 0.0, 4.24, -7752.82), new sun_con(13.0, -5.0, 2.01, 188491.91), new sun_con(13.0, 0.0, 2.65, 207.81), new sun_con(13.0, 0.0, 4.98, 29424.63), new sun_con(12.0, 0.0, 0.93, -7.99), new sun_con(10.0, 0.0, 2.21, 46941.14), new sun_con(10.0, 0.0, 3.59, -68.29), new sun_con(10.0, 0.0, 1.5, 21463.25), new sun_con(10.0, -9.0, 2.55, 157208.4) };

         ///*
         //Define the time unit 'u', measured in units of 10000 Julian years
         //from J2000.0.
         //*/

         u = System.Convert.ToDouble((jd - T0) / 3652500.0);

         ///*
         //Compute longitude and distance terms from the series.
         //*/

         for (i = 0; i <= 49; i++) {

            arg = System.Convert.ToDouble(con[i].alpha + con[i].nu * u);
            sum_lon += System.Convert.ToDouble(con[i].l * Sin(arg));
            sum_r += System.Convert.ToDouble(con[i].r * Cos(arg));
         }

         ///*
         //Compute longitude, latitude, and distance referred to mean equinox
         //and ecliptic of date.
         //*/

         lon = System.Convert.ToDouble(4.9353929 + 62833.196168 * u + factor * sum_lon);

         lon = lon % TWOPI;
         if (lon < 0.0) {
            lon += System.Convert.ToDouble(TWOPI);
         }

         lat = 0.0;

         dis = System.Convert.ToDouble(1.0001026 + factor * sum_r);

         ///*
         //Compute mean obliquity of the ecliptic.
         //*/

         t = u * 100.0;
         t2 = t * t;
         emean = System.Convert.ToDouble((0.001813 * t2 * t - 0.00059 * t2 - 46.815 * t + 84381.448) / RAD2SEC);

         ///*
         //Compute equatorial spherical coordinates referred to the mean equator
         //and equinox of date.
         //*/

         sin_lon = System.Convert.ToDouble(Sin(lon));
         ra = System.Convert.ToDouble(Atan2((Cos(emean) * sin_lon), Cos(lon)) * RAD2DEG);
         ra = ra % 360.0;
         if (ra < 0.0) {
            ra += 360.0;
         }
         ra = ra / 15.0;

         dec = System.Convert.ToDouble(Asin(Sin(emean) * sin_lon) * RAD2DEG);

      }

   }
}