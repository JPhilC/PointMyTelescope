using System;

namespace ASCOM.Astrometry
{
   sealed class DeltatCode
   {
      /// <summary>
      /// Calculates the value of DeltaT over a wide range of hstoric and future Julian dates
      /// </summary>
      /// <param name="tjd">Julian Date of interest</param>
      /// <returns>DelatT value at the given Julian date</returns>
      /// <remarks>
      /// Post 2011, calculation is effected throgh a 2nd order polynomial best fit to real DeltaT data from: http://maia.usno.navy.mil/ser7/deltat.data
      /// together with projections of DeltaT from: http://maia.usno.navy.mil/ser7/deltat.preds
      /// The analysis spreadsheets for DeltaT values at dates post 2011 are stored in the \NOVAS\DeltaT Predictions folder of the ASCOM source tree.
      /// </remarks>
      // VBConversions Note: Former VB static variables moved to class level because they aren't supported in C#.
      static double DeltaTCalc_ans = 0;
      static double DeltaTCalc_lasttjd = System.Convert.ToDouble(double.MinValue);

      static public double DeltaTCalc(double tjd)
      {

         const double TABSTART1620 = 1620.0;
         const int TABSIZ = 392;

         double YearFraction = 0;
         double B = 0;

         // static double ans = 0; VBConversions Note: Static variable moved to class level and renamed DeltaTCalc_ans. Local static variables are not supported in C#. // Variable To hold the last calculated DelatT value
         // static double lasttjd = System.Convert.ToDouble(double.MinValue); VBConversions Note: Static variable moved to class level and renamed DeltaTCalc_lasttjd. Local static variables are not supported in C#. // Initialise the last used Julian date to a value that will force a calculation on first time use

         // Performance optimisation:
         if (tjd == DeltaTCalc_lasttjd) {
            return DeltaTCalc_ans; // Return last calculated value if tjd is same as last call
         }
         DeltaTCalc_lasttjd = tjd; // Save the Julian date for use in the caching approach above

         YearFraction = System.Convert.ToDouble(2000.0 + (tjd - GlobalItems.T0) / 365.25); // This calculation is accurate enough for our purposes here (T0 = 2451545.0 is TDB Julian date of epoch J2000.0)

         // DATE RANGE January 2017 Onwards - The analysis was performed on 29th December 2016 and creates values within 0.12 of a second of the projections to Q3 2019
         if ((YearFraction >= 2017.0) && (YearFraction < double.MaxValue)) {
            DeltaTCalc_ans = System.Convert.ToDouble((0.02465436 * YearFraction * YearFraction) + (-98.92626556 * YearFraction) + 99301.85784308);
            return (DeltaTCalc_ans);
         }

         // DATE RANGE October 2015 Onwards - The analysis was performed on 24th October 2015 and creates values within 0.05 of a second of the projections to Q2 2018
         if ((YearFraction >= 2015.75) && (YearFraction < double.MaxValue)) {
            DeltaTCalc_ans = System.Convert.ToDouble((0.02002376 * YearFraction * YearFraction) + (-80.27921003 * YearFraction) + 80529.32);
            return (DeltaTCalc_ans);
         }

         // DATE RANGE October 2011 to September 2015 - The analysis was performed on 6th February 2014 and creates values within 0.2 of a second of the projections to Q1 2016
         if ((YearFraction >= 2011.75) && (YearFraction < 2015.75)) {
            DeltaTCalc_ans = System.Convert.ToDouble((0.00231189 * YearFraction * YearFraction) + (-8.85231952 * YearFraction) + 8518.54);
            return (DeltaTCalc_ans);
         }

         // DATE RANGE January 2011 to September 2011
         if ((YearFraction >= 2011.0) && (YearFraction < 2011.75)) {
            // Following now superseded by above for 2012-16, this is left in for consistency with previous behaviour
            // Use polynomial given at http://sunearth.gsfc.nasa.gov/eclipse/SEcat5/deltatpoly.html as retrtieved on 11-Jan-2009
            B = YearFraction - 2000.0;
            DeltaTCalc_ans = System.Convert.ToDouble(62.92 + (B * (0.32217 + (B * 0.005589))));
            return (DeltaTCalc_ans);
         }

         // Setup for pre 2011 calculations using Bob Denny's original code

         ///* Note, Stephenson and Morrison's table starts at the year 1630.
         // * The Chapronts' table does not agree with the Almanac prior to 1630.
         // * The actual accuracy decreases rapidly prior to 1780.
         // */
         //static short dt[] = {
         short[] dt = new short[] {(short)
            12400, (short) 11900, (short) 11500, (short) 11000, (short) 10600, (short) 10200, (short) 9800, (short) 9500, (short) 9100, (short) 8800, (short)
            8500, (short) 8200, (short) 7900, (short) 7700, (short) 7400, (short) 7200, (short) 7000, (short) 6700, (short) 6500, (short) 6300, (short)
            6200, (short) 6000, (short) 5800, (short) 5700, (short) 5500, (short) 5400, (short) 5300, (short) 5100, (short) 5000, (short) 4900, (short)
            4800, (short) 4700, (short) 4600, (short) 4500, (short) 4400, (short) 4300, (short) 4200, (short) 4100, (short) 4000, (short) 3800, (short)
            3700, (short) 3600, (short) 3500, (short) 3400, (short) 3300, (short) 3200, (short) 3100, (short) 3000, (short) 2800, (short) 2700, (short)
            2600, (short) 2500, (short) 2400, (short) 2300, (short) 2200, (short) 2100, (short) 2000, (short) 1900, (short) 1800, (short) 1700, (short)
            1600, (short) 1500, (short) 1400, (short) 1400, (short) 1300, (short) 1200, (short) 1200, (short) 1100, (short) 1100, (short) 1000, (short)
            1000, (short) 1000, (short) 900, (short) 900, (short) 900, (short) 900, (short) 900, (short) 900, (short) 900, (short) 900, (short)
            900, (short) 900, (short) 900, (short) 900, (short) 900, (short) 900, (short) 900, (short) 900, (short) 1000, (short) 1000, (short)
            1000, (short) 1000, (short) 1000, (short) 1000, (short) 1000, (short) 1000, (short) 1000, (short) 1100, (short) 1100, (short) 1100, (short)
            1100, (short) 1100, (short) 1100, (short) 1100, (short) 1100, (short) 1100, (short) 1100, (short) 1100, (short) 1100, (short) 1100, (short)
            1100, (short) 1100, (short) 1100, (short) 1100, (short) 1200, (short) 1200, (short) 1200, (short) 1200, (short) 1200, (short) 1200, (short)
            1200, (short) 1200, (short) 1200, (short) 1200, (short) 1300, (short) 1300, (short) 1300, (short) 1300, (short) 1300, (short) 1300, (short)
            1300, (short) 1400, (short) 1400, (short) 1400, (short) 1400, (short) 1400, (short) 1400, (short) 1400, (short) 1500, (short) 1500, (short)
            1500, (short) 1500, (short) 1500, (short) 1500, (short) 1500, (short) 1600, (short) 1600, (short) 1600, (short) 1600, (short) 1600, (short)
            1600, (short) 1600, (short) 1600, (short) 1600, (short) 1600, (short) 1700, (short) 1700, (short) 1700, (short) 1700, (short) 1700, (short)
            1700, (short) 1700, (short) 1700, (short) 1700, (short) 1700, (short) 1700, (short) 1700, (short) 1700, (short) 1700, (short) 1700, (short)
            1700, (short) 1700, (short) 1600, (short) 1600, (short) 1600, (short) 1600, (short) 1500, (short) 1500, (short) 1400, (short) 1400, (short)
            1370, (short) 1340, (short) 1310, (short) 1290, (short) 1270, (short) 1260, (short) 1250, (short) 1250, (short) 1250, (short) 1250, (short)
            1250, (short) 1250, (short) 1250, (short) 1250, (short) 1250, (short) 1250, (short) 1250, (short) 1240, (short) 1230, (short) 1220, (short)
            1200, (short) 1170, (short) 1140, (short) 1110, (short) 1060, (short) 1020, (short) 960, (short) 910, (short) 860, (short) 800, (short)
            750, (short) 700, (short) 660, (short) 630, (short) 600, (short) 580, (short) 570, (short) 560, (short) 560, (short) 560, (short)
            570, (short) 580, (short) 590, (short) 610, (short) 620, (short) 630, (short) 650, (short) 660, (short) 680, (short) 690, (short)
            710, (short) 720, (short) 730, (short) 740, (short) 750, (short) 760, (short) 770, (short) 770, (short) 780, (short) 780, (short)
            788, (short) 782, (short) 754, (short) 697, (short) 640, (short) 602, (short) 541, (short) 410, (short) 292, (short) 182, (short)
            161, (short) 10, (short) (-102), (short) (-128), (short) (-269), (short) (-324), (short) (-364), (short) (-454), (short) (-471), (short) (-511), (short)
            -540, (short) (-542), (short) (-520), (short) (-546), (short) (-546), (short) (-579), (short) (-563), (short) (-564), (short) (-580), (short) (-566), (short)
            -587, (short) (-601), (short) (-619), (short) (-664), (short) (-644), (short) (-647), (short) (-609), (short) (-576), (short) (-466), (short) (-374), (short)
            -272, (short) (-154), (short) (-2), (short) 124, (short) 264, (short) 386, (short) 537, (short) 614, (short) 775, (short) 913, (short)
            1046, (short) 1153, (short) 1336, (short) 1465, (short) 1601, (short) 1720, (short) 1824, (short) 1906, (short) 2025, (short) 2095, (short)
            2116, (short) 2225, (short) 2241, (short) 2303, (short) 2349, (short) 2362, (short) 2386, (short) 2449, (short) 2434, (short) 2408, (short)
            2402, (short) 2400, (short) 2387, (short) 2395, (short) 2386, (short) 2393, (short) 2373, (short) 2392, (short) 2396, (short) 2402, (short)
            2433, (short) 2483, (short) 2530, (short) 2570, (short) 2624, (short) 2677, (short) 2728, (short) 2778, (short) 2825, (short) 2871, (short)
            2915, (short) 2957, (short) 2997, (short) 3036, (short) 3072, (short) 3107, (short) 3135, (short) 3168, (short) 3218, (short) 3268, (short)
            3315, (short) 3359, (short) 3400, (short) 3447, (short) 3503, (short) 3573, (short) 3654, (short) 3743, (short) 3829, (short) 3920, (short)
            4018, (short) 4117, (short) 4223, (short) 4337, (short) 4449, (short) 4548, (short) 4646, (short) 4752, (short) 4853, (short) 4959, (short)
            5054, (short) 5138, (short) 5217, (short) 5296, (short) 5379, (short) 5434, (short) 5487, (short) 5532, (short) 5582, (short) 5630, (short)
            5686, (short) 5757, (short) 5831, (short) 5912, (short) 5998, (short) 6078, (short) 6163, (short) 6230, (short) 6296, (short) 6347, (short)
            6383, (short) 6409, (short) 6430, (short) 6447, (short) 6457, (short) 6469, (short) 6485, (short) 6515, (short) 6546, (short) 6570, (short)
            6650, (short) 6710};
         // Change TABEND and TABSIZ if you add/delete anything

         // Calculate  DeltaT = ET - UT in seconds.  Describes the irregularities of the Earth rotation rate in the ET time scale.
         double p = 0;
         int[] d = new int[7];
         int i = 0;
         int iy = 0;
         int k = 0;

         // DATE RANGE <1620
         if (YearFraction < TABSTART1620) {
            if (YearFraction >= 948.0) {
               ///* Stephenson and Morrison, stated domain is 948 to 1600:
               // * 25.5(centuries from 1800)^2 - 1.9159(centuries from 1955)^2
               // */
               B = System.Convert.ToDouble(0.01 * (YearFraction - 2000.0));
               DeltaTCalc_ans = System.Convert.ToDouble((23.58 * B + 100.3) * B + 101.6);
            }
            else {
               ///* Borkowski */
               B = System.Convert.ToDouble(0.01 * (YearFraction - 2000.0) + 3.75);
               DeltaTCalc_ans = System.Convert.ToDouble(35.0 * B * B + 40.0);
            }
            return DeltaTCalc_ans;
         }

         //DATE RANGE 1620 to 2011

         // Besselian interpolation from tabulated values. See AA page K11.
         // Index into the table.
         p = System.Convert.ToDouble(Math.Floor(YearFraction));
         iy = (int)(p - TABSTART1620); //// rbd - added cast
                                       ///* Zeroth order estimate is value at start of year */
         DeltaTCalc_ans = dt[iy];
         k = iy + 1;
         if (k >= TABSIZ) {
            goto done; // /* No data, can't go on. */
         }

         ///* The fraction of tabulation interval */
         p = YearFraction - p;

         ///* First order interpolated value */
         DeltaTCalc_ans += p * (dt[k] - dt[iy]);
         if ((iy - 1 < 0) || (iy + 2 >= TABSIZ)) {
            goto done; // /* can't do second differences */
         }

         ///* Make table of first differences */
         k = iy - 2;
         for (i = 0; i <= 4; i++) {
            if ((k < 0) || (k + 1 >= TABSIZ)) {
               d[i] = 0;
            }
            else {
               d[i] = dt[k + 1] - dt[k];
            }
            k++;
         }
         ///* Compute second differences */
         for (i = 0; i <= 3; i++) {
            d[i] = d[i + 1] - d[i];
         }
         B = System.Convert.ToDouble(0.25 * p * (p - 1.0));
         DeltaTCalc_ans += B * (d[1] + d[2]);
         if (iy + 2 >= TABSIZ) {
            goto done;
         }

         ///* Compute third differences */
         for (i = 0; i <= 2; i++) {
            d[i] = d[i + 1] - d[i];
         }
         B = System.Convert.ToDouble(2.0 * B / 3.0);
         DeltaTCalc_ans += System.Convert.ToDouble((p - 0.5) * B * d[1]);
         if ((iy - 2 < 0) || (iy + 3 > TABSIZ)) {
            goto done;
         }

         ///* Compute fourth differences */
         for (i = 0; i <= 1; i++) {
            d[i] = d[i + 1] - d[i];
         }
         B = System.Convert.ToDouble(0.125 * B * (p + 1.0) * (p - 2.0));
         DeltaTCalc_ans += B * (d[0] + d[1]);

         done:
         ///* Astronomical Almanac table is corrected by adding the expression
         // *     -0.000091 (ndot + 26)(year-1955)^2  seconds
         // * to entries prior to 1955 (AA page K8), where ndot is the secular
         // * tidal term in the mean motion of the Moon.
         // *
         // * Entries after 1955 are referred to atomic time standards and
         // * are not affected by errors in Lunar or planetary theory.
         // */
         DeltaTCalc_ans *= 0.01;
         if (YearFraction < 1955.0) {
            B = YearFraction - 1955.0;
            DeltaTCalc_ans += System.Convert.ToDouble(-0.000091 * (-25.8 + 26.0) * B * B);
         }
         return DeltaTCalc_ans;

      }
   }
}