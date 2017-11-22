using System.Runtime.InteropServices;
using ASCOM.Utilities.Exceptions;
using System;

namespace ASCOM.Astrometry.SOFA
{
   /// <summary>
   /// This class presents a subset of the SOFA (Standards of Fundamental Astronomy) astrometry routines in a form that is easily accessible from both 32bit and 64bit .NET and
   /// COM applications.
   /// </summary>
   /// <remarks>
   /// SOFA operates under the auspices of the International Astronomical Union (IAU) to provide algorithms and software for use in astronomical computing. The entire SOFA
   /// collection comprises many authoritative routines across a number of areas including:
   /// <list type="bullet">
   /// <item><description>Astrometry</description></item>
   /// <item><description>Calendars</description></item>
   /// <item><description>Time Scales</description></item>
   /// <item><description>Earth rotation and sidereal time</description></item>
   /// <item><description>Ephemerides (medium precision)</description></item>
   /// <item><description>Geocentric/geodetic transformations</description></item>
   /// <item><description>Precession, nutation, polar motion</description></item>
   /// <item><description>Star space motion</description></item>
   /// <item><description>Star catalogue conversion</description></item>
   /// </list>
   /// <para><b>The class's functionality is provided by underlying DLLs compiled from unmodified original C source code distributed by SOFA but the class does not constitute software provided by and/or endorsed by SOFA.
   /// No change whatsoever has been made to the algorithms implemented by SOFA that realise IAU standards.</b>
   /// </para>
   /// <para>SOFA provides a validation routine to confirm that the compiled library provides expected results. 32 and 64bit versions of this routine (SofaTestXX.exe and SofaTestXX-64.exe, where XX is the issue number)
   /// are included in this distribution and can be found in the Common Files\ASCOM\Astrometry directory. To run them open a command prompt in the Astrometry directory and enter the commands SofaTest10 /verbose and SofaTest10-64 /verbose.
   /// The susbset of these tests that is relevant to the routines presented in this component have also been incorporated in the ASCOM Diagnostics tool and expected operation of the SOFA routnines can be confirmed through this tool.</para>
   /// <para>Further information on the full library of SOFA routines is available here: http://www.iausofa.org/ </para>
   /// </remarks>
   public class SOFA : IDisposable
   {

      private const string SOFA32DLL = "SOFA12.dll"; //Names of SOFA 32 and 64bit DLL files
      private const string SOFA64DLL = "SOFA12-64.dll";

      // Release and revision constants
      private const int SOFA_RELEASE_NUMBER = 12;
      private const string SOFA_ISSUE_DATE = "2016-05-03";
      private const int SOFA_REVISION_NUMBER = 3; // Not presented in the interface, maintained here for reference
      private const string SOFA_REVISION_DATE = "2016-12-23";

      private IntPtr SofaDllHandle;

      #region New and IDisposable
      /// <summary>
      /// Creates a new instance of the SOFA component
      /// </summary>
      /// <exception cref="HelperException">Thrown if the SOFA support library DLL cannot be loaded</exception>
      /// <remarks></remarks>
      public SOFA()
      {
         bool rc;
         string SofaDllFile = "";
         System.Text.StringBuilder ReturnedPath = new System.Text.StringBuilder(260);
         int LastError = 0;



         //Find the root location of the common files directory containing the ASCOM support files.
         //On a 32bit system this is \Program Files\Common Files
         //On a 64bit system this is \Program Files (x86)\Common Files
         string exeLocation = System.Reflection.Assembly.GetExecutingAssembly().CodeBase;
         if (Is64Bit()) // 64bit application so find the 32bit folder location
         {
            rc = SHGetSpecialFolderPath(IntPtr.Zero, ReturnedPath, CSIDL_PROGRAM_FILES_COMMONX86, false);
            SofaDllFile = System.IO.Path.Combine(exeLocation, SOFA64DLL);
         }
         else //32bit application so just go with the .NET returned value
         {
            SofaDllFile = System.IO.Path.Combine(exeLocation, SOFA32DLL);
         }


         SofaDllHandle = LoadLibrary(SofaDllFile);
         LastError = Marshal.GetLastWin32Error();

         if (SofaDllHandle == IntPtr.Zero) // Loaded successfully
         {
            throw new DllNotFoundException("Error code returned from LoadLibrary when loading SOFA library: " + LastError.ToString("X8"));
         }
      }

      private bool disposedValue = false; // To detect redundant calls

      // IDisposable
      protected virtual void Dispose(bool disposing)
      {
         if (!this.disposedValue) {
            // Free your own state (unmanaged objects) and set large fields to null.
            try {
               FreeLibrary(SofaDllHandle);
            }
            catch {
            } // Free the SOFA library but don't return any error value

         }
         this.disposedValue = true;
      }

      // This code added by Visual Basic to correctly implement the disposable pattern.
      /// <summary>
      /// Cleans up the SOFA object
      /// </summary>
      /// <remarks></remarks>
      public void Dispose()
      {
         // Do not change this code.  Put cleanup code in Dispose(ByVal disposing As Boolean) above.
         Dispose(true);
         GC.SuppressFinalize(this);
      }
      #endregion

      #region Public SOFA Interface - SOFA Members
      /// <summary>
      /// Major number of the SOFA issue currently used by this component.
      /// </summary>
      /// <returns>Integer issue number</returns>
      /// <remarks></remarks>
      public int SofaReleaseNumber()
      {
         return SOFA_RELEASE_NUMBER; // The release number of the issue being used
      }

      /// <summary>
      /// Release date of the SOFA issue currently used by this component.
      /// </summary>
      /// <returns>String date in format yyyy-mm-dd</returns>
      /// <remarks></remarks>
      public string SofaIssueDate()
      {
         return SOFA_ISSUE_DATE; // Release date of the fundamental software issue currently being used
      }

      /// <summary>
      /// Release date of the revision to the SOFA Issue that is actually being used by this component.
      /// </summary>
      /// <returns>String date in format yyyy-mm-dd</returns>
      /// <remarks>When a new issue is employed that doesn't yet have a revision, this mehtod will return the SofaIssueDate</remarks>
      public string SofaRevisionDate()
      {
         return SOFA_REVISION_DATE; // Release date of the revision currently being used
      }

      /// <summary>
      /// Convert degrees, arcminutes, arcseconds to radians.
      /// </summary>
      /// <param name="s">Sign:  '-' = negative, otherwise positive</param>
      /// <param name="ideg">Degrees</param>
      /// <param name="iamin">Arcminutes</param>
      /// <param name="asec">Arcseconds</param>
      /// <param name="rad">Angle in radian</param>
      /// <returns>Status:  0 = OK, 1 = ideg outside range 0-359, 2 = iamin outside range 0-59, 3 = asec outside range 0-59.999...</returns>
      /// <remarks>
      /// Notes:
      /// <list type="number">
      /// <item><description>The result is computed even if any of the range checks fail.</description></item>
      /// <item><description>Negative ideg, iamin and/or asec produce a warning status, but the absolute value is used in the conversion.</description></item>
      /// <item><description>If there are multiple errors, the status value reflects only the first, the smallest taking precedence.</description></item>
      /// </list>
      /// </remarks>
      public int Af2a(string s, int ideg, int iamin, double asec, ref double rad)
      {

         short RetCode = 0;

         if (string.IsNullOrEmpty(s)) {
            s = " "; // Fix any invalid sign values
         }

         if (Is64Bit()) {
            Af2a64(s.ToCharArray()[0], System.Convert.ToInt16(Convert.ToInt16(ideg)), System.Convert.ToInt16(Convert.ToInt16(iamin)), asec, ref rad);
         }
         else {
            Af2a32(s.ToCharArray()[0], System.Convert.ToInt16(Convert.ToInt16(ideg)), System.Convert.ToInt16(Convert.ToInt16(iamin)), asec, ref rad);
         }

         return Convert.ToInt32(RetCode);
      }

      /// <summary>
      /// Normalize angle into the range 0 &lt;= a &lt; 2pi.
      /// </summary>
      /// <param name="a">Angle (radians)</param>
      /// <returns>Angle in range 0-2pi</returns>
      /// <remarks></remarks>
      public double Anp(double a)
      {
         double RetVal = 0;

         if (Is64Bit()) {
            RetVal = Anp64(a);
         }
         else {
            RetVal = Anp32(a);
         }

         return RetVal;
      }

      /// <summary>
      /// Transform ICRS star data, epoch J2000.0, to CIRS using the SOFA Atci13 function.
      /// </summary>
      /// <param name="rc">ICRS right ascension at J2000.0 (radians, Note 1)</param>
      /// <param name="dc">ICRS declination at J2000.0 (radians, Note 1)</param>
      /// <param name="pr">RA proper motion (radians/year; Note 2)</param>
      /// <param name="pd">Dec proper motion (radians/year)</param>
      /// <param name="px">parallax (arcsec)</param>
      /// <param name="rv">radial velocity (km/s, +ve if receding)</param>
      /// <param name="date1">TDB as a 2-part Julian Date (Note 3)</param>
      /// <param name="date2">TDB as a 2-part Julian Date (Note 3)</param>
      /// <param name="ri">CIRS geocentric RA (radians)</param>
      /// <param name="di">CIRS geocentric Dec (radians)</param>
      /// <param name="eo">equation of the origins (ERA-GST, Note 5)</param>
      /// <remarks>
      /// Notes:
      /// <list type="number">
      /// <item><description>Star data for an epoch other than J2000.0 (for example from the Hipparcos catalog, which has an epoch of J1991.25) will require a preliminary call to iauPmsafe before use.</description></item>
      /// <item><description>The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.</description></item>
      /// <item><description> The TDB date date1+date2 is a Julian Date, apportioned in any convenient way between the two arguments.  For example, JD(TDB)=2450123.8g could be expressed in any of these ways, among others:
      /// <table style="width:340px;" cellspacing="0">
      /// <col style="width:80px;"></col>
      /// <col style="width:80px;"></col>
      /// <col style="width:180px;"></col>
      /// <tr>
      /// <td colspan="1" align="center" rowspan="1" style="width: 80px; padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-top-color: #000000; border-top-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px;
      /// background-color: #00ffff;" width="110px">
      /// <b>Date 1</b></td>
      /// <td colspan="1" rowspan="1" align="center" style="width: 80px; padding-right: 10px; padding-left: 10px;
      /// border-top-color: #000000; border-top-style: Solid;
      /// border-right-style: Solid; border-right-color: #000000;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px;
      /// background-color: #00ffff;" width="110px">
      /// <b>Date 2</b></td>
      /// <td colspan="1" rowspan="1" align="center" style="width: 180px; padding-right: 10px; padding-left: 10px;
      /// border-top-color: #000000; border-top-style: Solid;
      /// border-right-style: Solid; border-right-color: #000000;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px;
      /// background-color: #00ffff;" width="220px">
      /// <b>Method</b></td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      ///  2450123.8</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 0.0</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// JD method</td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 2451545.0</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// -1421.3</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// J2000 method</td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 2400000.5</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 50123.2</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// MJD method</td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 2450123.5</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 0.2</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// Date and time method</td>
      /// </tr>
      /// </table>
      /// <para>The JD method is the most natural and convenient to use in cases where the loss of several decimal digits of resolution is acceptable.  The J2000 method is best matched to the way the argument is handled internally
      /// and will deliver the optimum resolution.  The MJD method and the date and time methods are both good compromises between resolution and convenience.  For most applications of this function the choice will not be at all critical.</para>
      /// <para>TT can be used instead of TDB without any significant impact on accuracy.</para>
      /// </description></item>
      /// <item><description>The available accuracy is better than 1 milliarcsecond, limited mainly by the precession-nutation model that is used, namely IAU 2000A/2006.  Very close to solar system bodies, additional
      /// errors of up to several milliarcseconds can occur because of unmodeled light deflection;  however, the Sun's contribution is taken into account, to first order.  The accuracy limitations of
      /// the SOFA function iauEpv00 (used to compute Earth position and velocity) can contribute aberration errors of up to 5 microarcseconds.  Light deflection at the Sun's limb is uncertain at the 0.4 mas level.</description></item>
      /// <item><description>Should the transformation to (equinox based) apparent place be required rather than (CIO based) intermediate place, subtract the equation of the origins from the returned right ascension:
      /// RA = RI - EO. (The Anp function can then be applied, as required, to keep the result in the conventional 0-2pi range.)</description></item>
      /// </list>
      /// </remarks>
      public void Atci13(double rc, double dc, double pr, double pd, double px, double rv, double date1, double date2, ref double ri, ref double di, ref double eo)
      {
         this.CelestialToIntermediate(rc, dc, pr, pd, px, rv, date1, date2, ref ri, ref di, ref eo);
      }

      public void CelestialToIntermediate(double rc, double dc, double pr, double pd, double px, double rv, double date1, double date2, ref double ri, ref double di, ref double eo)
      {

         if (Is64Bit()) {
            Atci1364(rc, dc, pr, pd, px, rv, date1, date2, ref ri, ref di, ref eo);
         }
         else {
            Atci1332(rc, dc, pr, pd, px, rv, date1, date2, ref ri, ref di, ref eo);
         }
      }

      /// <summary>
      /// ICRS RA,Dec to observed place using the SOFA Atco13 function.
      /// </summary>
      /// <param name="rc">ICRS RA (radians, note 1)</param>
      /// <param name="dc">ICRS Dec (radians, note 2)</param>
      /// <param name="pr">RA Proper motion (radians/year)</param>
      /// <param name="pd">Dec Proper motion (radians/year</param>
      /// <param name="px">Parallax (arcsec)</param>
      /// <param name="rv">Radial veolcity (Km/s, +ve if receding</param>
      /// <param name="utc1">UTC Julian date (part 1, notes 3,4)</param>
      /// <param name="utc2">UTC Julian date (part 2, notes 3,4)</param>
      /// <param name="dut1">UT1 - UTC (seonds, note 5)</param>
      /// <param name="elong">Site longitude (radians, note 6)</param>
      /// <param name="phi">Site Latitude (radians, note 6)</param>
      /// <param name="hm">Site Height (meters, notes 6,8)</param>
      /// <param name="xp">Polar motion co-ordinate (radians, note 7)</param>
      /// <param name="yp">Polar motion co-ordinate (radians,note 7)</param>
      /// <param name="phpa">Site Presure (hPa = mB, note 8)</param>
      /// <param name="tc">Site Temperature (C)</param>
      /// <param name="rh">Site relative humidity (fraction in the range: 0.0 to 1.0)</param>
      /// <param name="wl">Observation wavem=length (micrometres, note 9)</param>
      /// <param name="aob">Observed Azimuth (radians)</param>
      /// <param name="zob">Observed Zenith distance (radians)</param>
      /// <param name="hob">Observed Hour Angle (radians)</param>
      /// <param name="dob">Observed Declination (radians)</param>
      /// <param name="rob">Observed RA (radians)</param>
      /// <param name="eo">Equation of the origins (ERA-GST)</param>
      /// <returns>+1 = dubious year (Note 4), 0 = OK, -1 = unacceptable date</returns>
      /// <remarks>
      /// Notes:
      /// <list type="number">
      /// <item><description>Star data for an epoch other than J2000.0 (for example from the Hipparcos catalog, which has an epoch of J1991.25) will require a preliminary call to iauPmsafe before use.</description></item>
      /// <item><description>The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.</description></item>
      /// <item><description>utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any convenient way between the two arguments, for example where utc1 is the Julian Day Number and utc2 is the fraction of a day.
      /// <para>However, JD cannot unambiguously represent UTC during a leap second unless special measures are taken.  The convention in the present function is that the JD day represents UTC days whether the length is 86399, 86400 or 86401 SI seconds.</para>
      /// <para>Applications should use the function iauDtf2d to convert from calendar date and time of day into 2-part quasi Julian Date, as it implements the leap-second-ambiguity convention just described.</para></description></item>
      /// <item><description>The warning status "dubious year" flags UTCs that predate the introduction of the time scale or that are too far in the future to be trusted.  See iauDat for further details.</description></item>
      /// <item><description>UT1-UTC is tabulated in IERS bulletins.  It increases by exactly one second at the end of each positive UTC leap second, introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This practice is under review, and in the future UT1-UTC may grow essentially without limit.</description></item>
      /// <item><description>The geographical coordinates are with respect to the WGS84 reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the longitude required by the present function is east-positive (i.e. right-handed), in accordance with geographical convention.</description></item>
      /// <item><description>The polar motion xp,yp can be obtained from IERS bulletins.  The values are the coordinates (in radians) of the Celestial Intermediate Pole with respect to the International Terrestrial Reference System (see IERS Conventions 2003), measured along the meridians 0 and 90 deg west respectively.  For many applications, xp and yp can be set to zero.</description></item>
      /// <item><description>If hm, the height above the ellipsoid of the observing station in meters, is not known but phpa, the pressure in hPa (=mB), is available, an adequate estimate of hm can be obtained from the expression:
      /// <p style="margin-left:25px;font-family:Lucida Conosle,Monospace"><b>hm = -29.3 * tsl * log ( phpa / 1013.25 );</b></p>
      /// <para>where tsl is the approximate sea-level air temperature in K (See Astrophysical Quantities, C.W.Allen, 3rd edition, section 52).  Similarly, if the pressure phpa is not known, it can be estimated from the height of the observing station, hm, as follows:</para>
      /// <p style="margin-left:25px;font-family:Lucida Conosle,Monospace"><b>phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );</b></p>
      /// <para>Note, however, that the refraction is nearly proportional to the pressure and that an accurate phpa value is important for precise work.</para></description></item>
      /// <item><description>The argument wl specifies the observing wavelength in micrometers.  The transition from optical to radio is assumed to occur at 100 micrometers (about 3000 GHz).</description></item>
      /// <item><description>The accuracy of the result is limited by the corrections for refraction, which use a simple A*tan(z) + B*tan^3(z) model. Providing the meteorological parameters are known accurately and there are no gross local effects, the predicted observed coordinates should be within 0.05 arcsec (optical) or 1 arcsec (radio) for a zenith distance of less than 70 degrees, better than 30 arcsec (optical or radio) at 85 degrees and better than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.
      /// <para>Without refraction, the complementary functions iauAtco13 and iauAtoc13 are self-consistent to better than 1 microarcsecond all over the celestial sphere.  With refraction included, consistency falls off at high zenith distances, but is still better than 0.05 arcsec at 85 degrees.</para></description></item>
      /// <item><description>"Observed" Az,ZD means the position that would be seen by a perfect geodetically aligned theodolite.  (Zenith distance is used rather than altitude in order to reflect the fact that no allowance is made for depression of the horizon.)  This is related to the observed HA,Dec via the standard rotation, using the geodetic latitude (corrected for polar motion), while the observed HA and RA are related simply through the Earth rotation angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus means the position that would be seen by a perfect equatorial with its polar axis aligned to the Earth's axis of rotation.</description></item>
      /// <item><description>It is advisable to take great care with units, as even unlikely values of the input parameters are accepted and processed in accordance with the models used.</description></item>
      /// </list>
      /// </remarks>
      public int Atco13(double rc, double dc, double pr, double pd, double px, double rv, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double aob, ref double zob, ref double hob, ref double dob, ref double rob, ref double eo)
      {
         return this.CelestialToObserved(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref aob, ref zob, ref hob, ref dob, ref rob, ref eo);
      }

      public int CelestialToObserved(double rc, double dc, double pr, double pd, double px, double rv, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double aob, ref double zob, ref double hob, ref double dob, ref double rob, ref double eo)
      {

         short RetCode = 0;

         if (Is64Bit()) {
            RetCode = Atco1364(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref aob, ref zob, ref hob, ref dob, ref rob, ref eo);
         }
         else {
            RetCode = Atco1332(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref aob, ref zob, ref hob, ref dob, ref rob, ref eo);
         }

         return Convert.ToInt32(RetCode);
      }

      /// <summary>
      /// Encode date and time fields into 2-part Julian Date (or in the case of UTC a quasi-JD form that includes special provision for leap seconds).
      /// </summary>
      /// <param name="scale">Time scale ID (Note 1)</param>
      /// <param name="iy">Year in Gregorian calendar (Note 2)</param>
      /// <param name="im">Month in Gregorian calendar (Note 2)</param>
      /// <param name="id">Day in Gregorian calendar (Note 2)</param>
      /// <param name="ihr">Hour</param>
      /// <param name="imn">Minute</param>
      /// <param name="sec">Seconds</param>
      /// <param name="d1">2-part Julian Date (Notes 3, 4)</param>
      /// <param name="d2">2-part Julian Date (Notes 3, 4)</param>
      /// <returns>Status: +3 = both of next two, +2 = time is after end of day (Note 5), +1 = dubious year (Note 6), 0 = OK, -1 = bad year, -2 = bad month, -3 = bad day, -4 = bad hour, -5 = bad minute, -6 = bad second (&lt;0)</returns>
      /// <remarks>
      /// Notes:
      /// <list type="number">
      /// <item><description>Scale identifies the time scale.  Only the value "UTC" (in upper case) is significant, and enables handling of leap seconds (see Note 4).</description></item>
      /// <item><description>For calendar conventions and limitations, see iauCal2jd.</description></item>
      /// <item><description>The sum of the results, d1+d2, is Julian Date, where normally d1 is the Julian Day Number and d2 is the fraction of a day.  In the case of UTC, where the use of JD is problematical, special conventions apply:  see the next note.</description></item>
      /// <item><description>JD cannot unambiguously represent UTC during a leap second unless special measures are taken.  The SOFA internal convention is that the quasi-JD day represents UTC days whether the length is 86399,
      /// 86400 or 86401 SI seconds.  In the 1960-1972 era there were smaller jumps (in either direction) each time the linear UTC(TAI) expression was changed, and these "mini-leaps" are also included in the SOFA convention.</description></item>
      /// <item><description>The warning status "time is after end of day" usually means that the sec argument is greater than 60.0.  However, in a day ending in a leap second the limit changes to 61.0 (or 59.0 in the case of a negative leap second).</description></item>
      /// <item><description>The warning status "dubious year" flags UTCs that predate the introduction of the time scale or that are too far in the future to be trusted.  See iauDat for further details.</description></item>
      /// <item><description>Only in the case of continuous and regular time scales (TAI, TT, TCG, TCB and TDB) is the result d1+d2 a Julian Date, strictly speaking.  In the other cases (UT1 and UTC) the result must be
      /// used with circumspection;  in particular the difference between two such results cannot be interpreted as a precise time interval.</description></item>
      /// </list>
      /// </remarks>
      public int Dtf2d(string scale, int iy, int im, int id, int ihr, int imn, double sec, ref double d1, ref double d2)
      {

         short RetCode = 0;

         if (Is64Bit()) {
            RetCode = Dtf2d64(scale, iy, im, id, ihr, imn, sec, ref d1, ref d2);
         }
         else {
            RetCode = Dtf2d32(scale, iy, im, id, ihr, imn, sec, ref d1, ref d2);
         }

         return Convert.ToInt32(RetCode);
      }

      /// <summary>
      /// Equation of the origins, IAU 2006 precession and IAU 2000A nutation.
      /// </summary>
      /// <param name="date1">TT as a 2-part Julian Date (Note 1)</param>
      /// <param name="date2">TT as a 2-part Julian Date (Note 1)</param>
      /// <returns>Equation of the origins in radians (Note 2)</returns>
      /// <remarks>
      /// Notes:
      /// <list type="number">
      /// <item><description> The TT date date1+date2 is a Julian Date, apportioned in any convenient way between the two arguments.  For example, JD(TT)=2450123.7 could be expressed in any of these ways, among others:
      /// <table style="width:340px;" cellspacing="0">
      /// <col style="width:80px;"></col>
      /// <col style="width:80px;"></col>
      /// <col style="width:180px;"></col>
      /// <tr>
      /// <td colspan="1" align="center" rowspan="1" style="width: 80px; padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-top-color: #000000; border-top-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px;
      /// background-color: #00ffff;" width="110px">
      /// <b>Date 1</b></td>
      /// <td colspan="1" rowspan="1" align="center" style="width: 80px; padding-right: 10px; padding-left: 10px;
      /// border-top-color: #000000; border-top-style: Solid;
      /// border-right-style: Solid; border-right-color: #000000;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px;
      /// background-color: #00ffff;" width="110px">
      /// <b>Date 2</b></td>
      /// <td colspan="1" rowspan="1" align="center" style="width: 180px; padding-right: 10px; padding-left: 10px;
      /// border-top-color: #000000; border-top-style: Solid;
      /// border-right-style: Solid; border-right-color: #000000;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px;
      /// background-color: #00ffff;" width="220px">
      /// <b>Method</b></td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      ///  2450123.7</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 0.0</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// JD method</td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 2451545.0</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// -1421.3</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// J2000 method</td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 2400000.5</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 50123.2</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// MJD method</td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 2450123.5</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 0.2</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// Date and time method</td>
      /// </tr>
      /// </table>
      /// <para>The JD method is the most natural and convenient to use in cases where the loss of several decimal digits of resolution is acceptable.  The J2000 method is best matched to the way the argument is handled internally
      /// and will deliver the optimum resolution.  The MJD method and the date and time methods are both good compromises between resolution and convenience.  For most applications of this function the choice will not be at all critical.</para>
      /// </description></item>
      /// <item><description> The equation of the origins is the distance between the true equinox and the celestial intermediate origin and, equivalently, the difference between Earth rotation angle and Greenwich
      /// apparent sidereal time (ERA-GST).  It comprises the precession (since J2000.0) in right ascension plus the equation of the equinoxes (including the small correction terms).</description></item>
      /// </list>
      /// </remarks>
      public double Eo06a(double date1, double date2)
      {

         double RetVal = 0;

         if (Is64Bit()) {
            RetVal = Eo06a64(date1, date2);
         }
         else {
            RetVal = Eo06a32(date1, date2);
         }

         return RetVal;
      }

      /// <summary>
      /// Transform star RA,Dec from geocentric CIRS to ICRS astrometric using the SOFA Atic13 function.
      /// </summary>
      /// <param name="ri">CIRS geocentric RA (radians)</param>
      /// <param name="di">CIRS geocentric Dec (radians)</param>
      /// <param name="date1">TDB as a 2-part Julian Date (Note 1)</param>
      /// <param name="date2">TDB as a 2-part Julian Date (Note 1)</param>
      /// <param name="rc">ICRS astrometric RA (radians)</param>
      /// <param name="dc">ICRS astrometric Dec (radians)</param>
      /// <param name="eo">equation of the origins (ERA-GST, Note 4)</param>
      /// <remarks>
      /// Notes:
      /// <list type="number">
      /// <item><description> The TDB date date1+date2 is a Julian Date, apportioned in any convenient way between the two arguments.  For example, JD(TDB)=2450123.8g could be expressed in any of these ways, among others:
      /// <table style="width:340px;" cellspacing="0">
      /// <col style="width:80px;"></col>
      /// <col style="width:80px;"></col>
      /// <col style="width:180px;"></col>
      /// <tr>
      /// <td colspan="1" align="center" rowspan="1" style="width: 80px; padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-top-color: #000000; border-top-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px;
      /// background-color: #00ffff;" width="110px">
      /// <b>Date 1</b></td>
      /// <td colspan="1" rowspan="1" align="center" style="width: 80px; padding-right: 10px; padding-left: 10px;
      /// border-top-color: #000000; border-top-style: Solid;
      /// border-right-style: Solid; border-right-color: #000000;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px;
      /// background-color: #00ffff;" width="110px">
      /// <b>Date 2</b></td>
      /// <td colspan="1" rowspan="1" align="center" style="width: 180px; padding-right: 10px; padding-left: 10px;
      /// border-top-color: #000000; border-top-style: Solid;
      /// border-right-style: Solid; border-right-color: #000000;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px;
      /// background-color: #00ffff;" width="220px">
      /// <b>Method</b></td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      ///  2450123.8</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 0.0</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// JD method</td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 2451545.0</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// -1421.3</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// J2000 method</td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 2400000.5</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 50123.2</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// MJD method</td>
      /// </tr>
      /// <tr>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-left-color: #000000; border-left-style: Solid;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 2450123.5</td>
      /// <td align="right" style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// 0.2</td>
      /// <td style="padding-right: 10px; padding-left: 10px;
      /// border-right-color: #000000; border-right-style: Solid;
      /// border-bottom-color: #000000; border-bottom-style: Solid;
      /// border-right-width: 1px; border-left-width: 1px; border-top-width: 1px; border-bottom-width: 1px; ">
      /// Date and time method</td>
      /// </tr>
      /// </table>
      /// <para>The JD method is the most natural and convenient to use in cases where the loss of several decimal digits of resolution is acceptable.  The J2000 method is best matched to the way the argument is handled internally
      /// and will deliver the optimum resolution.  The MJD method and the date and time methods are both good compromises between resolution and convenience.  For most applications of this function the choice will not be at all critical.</para>
      /// <para>TT can be used instead of TDB without any significant impact on accuracy.</para>
      /// </description></item>
      /// <item><description>Iterative techniques are used for the aberration and light deflection corrections so that the functions Atic13 and Atci13 are accurate inverses;
      /// even at the edge of the Sun's disk the discrepancy is only about 1 nanoarcsecond.</description></item>
      /// <item><description>The available accuracy is better than 1 milliarcsecond, limited mainly by the precession-nutation model that is used, namely IAU 2000A/2006.  Very close to solar system bodies, additional
      /// errors of up to several milliarcseconds can occur because of unmodeled light deflection;  however, the Sun's contribution is taken into account, to first order.  The accuracy limitations of
      /// the SOFA function iauEpv00 (used to compute Earth position and velocity) can contribute aberration errors of up to 5 microarcseconds.  Light deflection at the Sun's limb is uncertain at the 0.4 mas level.</description></item>
      /// <item><description>Should the transformation to (equinox based) J2000.0 mean place be required rather than (CIO based) ICRS coordinates, subtract the equation of the origins from the returned right ascension:
      /// RA = RI - EO.  (The Anp function can then be applied, as required, to keep the result in the conventional 0-2pi range.)</description></item>
      /// </list>
      /// </remarks>
      public void Atic13(double ri, double di, double date1, double date2, ref double rc, ref double dc, ref double eo)
      {
         this.IntermediateToCelestial(ri, di, date1, date2, ref rc, ref dc, ref eo);
      }

      public void IntermediateToCelestial(double ri, double di, double date1, double date2, ref double rc, ref double dc, ref double eo)
      {

         if (Is64Bit()) {
            Atic1364(ri, di, date1, date2, ref rc, ref dc, ref eo);
         }
         else {
            Atic1332(ri, di, date1, date2, ref rc, ref dc, ref eo);
         }
      }

      /// <summary>
      /// CIRS RA,Dec to observed place using the SOFA Atio13 funciton.
      /// </summary>
      /// <param name="ri">CIRS right ascension (CIO-based, radians)</param>
      /// <param name="di">CIRS declination (radians)</param>
      /// <param name="utc1">UTC as a 2-part quasi Julian Date (Notes 1,2)</param>
      /// <param name="utc2">UTC as a 2-part quasi Julian Date (Notes 1,2)</param>
      /// <param name="dut1">UT1-UTC (seconds, Note 3)</param>
      /// <param name="elong">longitude (radians, east +ve, Note 4)</param>
      /// <param name="phi">geodetic latitude (radians, Note 4)</param>
      /// <param name="hm">height above ellipsoid (m, geodetic Notes 4,6)</param>
      /// <param name="xp">polar motion coordinates (radians, Note 5)</param>
      /// <param name="yp">polar motion coordinates (radians, Note 5)</param>
      /// <param name="phpa">pressure at the observer (hPa = mB, Note 6)</param>
      /// <param name="tc">ambient temperature at the observer (deg C)</param>
      /// <param name="rh">relative humidity at the observer (range 0-1)</param>
      /// <param name="wl">wavelength (micrometers, Note 7)</param>
      /// <param name="aob">observed azimuth (radians: N=0,E=90)</param>
      /// <param name="zob">observed zenith distance (radians)</param>
      /// <param name="hob">observed hour angle (radians)</param>
      /// <param name="dob">observed declination (radians)</param>
      /// <param name="rob">observed right ascension (CIO-based, radians)</param>
      /// <returns> Status: +1 = dubious year (Note 2), 0 = OK, -1 = unacceptable date</returns>
      /// <remarks>
      /// <para>Notes:</para>
      /// <list type="number">
      /// <item><description>utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any convenient way between the two arguments, for example where utc1 is the Julian Day Number and utc2 is the fraction of a day.
      /// <para>However, JD cannot unambiguously represent UTC during a leap second unless special measures are taken.  The convention in the present function is that the JD day represents UTC days whether the length is 86399, 86400 or 86401 SI seconds.</para>
      /// <para>Applications should use the function iauDtf2d to convert from calendar date and time of day into 2-part quasi Julian Date, as it implements the leap-second-ambiguity convention just described.</para></description></item>
      /// <item><description>The warning status "dubious year" flags UTCs that predate the introduction of the time scale or that are too far in the future to be trusted.  See iauDat for further details.</description></item>
      /// <item><description>UT1-UTC is tabulated in IERS bulletins.  It increases by exactly one second at the end of each positive UTC leap second, introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This practice is under review, and in the future UT1-UTC may grow essentially without limit.</description></item>
      /// <item><description>The geographical coordinates are with respect to the WGS84 reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the longitude required by the present function is east-positive (i.e. right-handed), in accordance with geographical convention.</description></item>
      /// <item><description>The polar motion xp,yp can be obtained from IERS bulletins.  The values are the coordinates (in radians) of the Celestial Intermediate Pole with respect to the International Terrestrial
      /// Reference System (see IERS Conventions 2003), measured along the meridians 0 and 90 deg west respectively.  For many applications, xp and yp can be set to zero.</description></item>
      /// <item><description>If hm, the height above the ellipsoid of the observing station in meters, is not known but phpa, the pressure in hPa (=mB), is available, an adequate estimate of hm can be obtained from the expression:
      /// <p style="margin-left:25px;font-family:Lucida Conosle,Monospace"><b>hm = -29.3 * tsl * log ( phpa / 1013.25 );</b></p>
      /// <para>where tsl is the approximate sea-level air temperature in K (See Astrophysical Quantities, C.W.Allen, 3rd edition, section 52).  Similarly, if the pressure phpa is not known, it can be estimated from the height of the observing station, hm, as follows:</para>
      /// <p style="margin-left:25px;font-family:Lucida Conosle,Monospace"><b>phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );</b></p>
      /// <para>Note, however, that the refraction is nearly proportional to the pressure and that an accurate phpa value is important for precise work.</para></description></item>
      /// <item><description>The argument wl specifies the observing wavelength in micrometers.  The transition from optical to radio is assumed to occur at 100 micrometers (about 3000 GHz).</description></item>
      /// <item><description>"Observed" Az,ZD means the position that would be seen by a perfect geodetically aligned theodolite.  (Zenith distance is used rather than altitude in order to reflect the fact that no
      /// allowance is made for depression of the horizon.)  This is related to the observed HA,Dec via the standard rotation, using the geodetic latitude (corrected for polar motion), while the observed HA and RA are related simply through the Earth rotation
      /// angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus means the position that would be seen by a perfect equatorial with its polar axis aligned to the Earth's axis of rotation.</description></item>
      /// <item><description>The accuracy of the result is limited by the corrections for refraction, which use a simple A*tan(z) + B*tan^3(z) model. Providing the meteorological parameters are known accurately and there are no gross local effects, the predicted astrometric
      /// coordinates should be within 0.05 arcsec (optical) or 1 arcsec (radio) for a zenith distance of less than 70 degrees, better than 30 arcsec (optical or radio) at 85 degrees and better than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.</description></item>
      /// <item><description>The complementary functions iauAtio13 and iauAtoi13 are self-consistent to better than 1 microarcsecond all over the celestial sphere.</description></item>
      /// <item><description>It is advisable to take great care with units, as even unlikely values of the input parameters are accepted and processed in accordance with the models used.</description></item>
      /// </list>
      /// </remarks>
      public int Atio13(double ri, double di, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double aob, ref double zob, ref double hob, ref double dob, ref double rob)
      {
         return this.IntermediateToObserved(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref aob, ref zob, ref hob, ref dob, ref rob);
      }

      public int IntermediateToObserved(double ri, double di, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double aob, ref double zob, ref double hob, ref double dob, ref double rob)
      {
         short RetCode = 0;

         if (Is64Bit()) {
            RetCode = Atio1364(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref aob, ref zob, ref hob, ref dob, ref rob);
         }
         else {
            RetCode = Atio1332(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref aob, ref zob, ref hob, ref dob, ref rob);
         }

         return Convert.ToInt32(RetCode);
      }

      /// <summary>
      /// Observed place at a groundbased site to to ICRS astrometric RA,Dec using the SOFA Atoc13 function.
      /// </summary>
      /// <param name="type">type of coordinates - "R", "H" or "A" (Notes 1,2)</param>
      /// <param name="ob1">observed Az, HA or RA (radians; Az is N=0,E=90)</param>
      /// <param name="ob2"> observed ZD or Dec (radians)</param>
      /// <param name="utc1">UTC as a 2-part quasi Julian Date (Notes 3,4)</param>
      /// <param name="utc2">UTC as a 2-part quasi Julian Date (Notes 3,4)</param>
      /// <param name="dut1">UT1-UTC (seconds, Note 5)</param>
      /// <param name="elong">longitude (radians, east +ve, Note 6)</param>
      /// <param name="phi">geodetic latitude (radians, Note 6)</param>
      /// <param name="hm">height above ellipsoid (m, geodetic Notes 6,8)</param>
      /// <param name="xp">polar motion coordinates (radians, Note 7)</param>
      /// <param name="yp">polar motion coordinates (radians, Note 7)</param>
      /// <param name="phpa">pressure at the observer (hPa = mB, Note 8)</param>
      /// <param name="tc">ambient temperature at the observer (deg C)</param>
      /// <param name="rh">relative humidity at the observer (range 0-1)</param>
      /// <param name="wl">wavelength (micrometers, Note 9)</param>
      /// <param name="rc">ICRS astrometric RA (radians)</param>
      /// <param name="dc">ICRS astrometric Dec (radians)</param>
      /// <returns>Status: +1 = dubious year (Note 4), 0 = OK, -1 = unacceptable date</returns>
      /// <remarks>
      /// <para>Notes:</para>
      /// <list type="number">
      /// <item><description>"Observed" Az,ZD means the position that would be seen by a perfect geodetically aligned theodolite.  (Zenith distance is used rather than altitude in order to reflect the fact that no
      /// allowance is made for depression of the horizon.)  This is related to the observed HA,Dec via the standard rotation, using the geodetic latitude (corrected for polar motion), while the
      /// observed HA and RA are related simply through the Earth rotation angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus means the position that would be seen by a perfect equatorial with its polar axis aligned to the Earth's axis of rotation.</description></item>
      /// <item><description>Only the first character of the type argument is significant. "R" or "r" indicates that ob1 and ob2 are the observed right ascension and declination;  "H" or "h" indicates that they are hour angle (west +ve) and declination;  anything else ("A" or
      /// "a" is recommended) indicates that ob1 and ob2 are azimuth (north zero, east 90 deg) and zenith distance.</description></item>
      /// <item><description>utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any convenient way between the two arguments, for example where utc1 is the Julian Day Number and utc2 is the fraction of a day.
      /// <para>However, JD cannot unambiguously represent UTC during a leap second unless special measures are taken.  The convention in the present function is that the JD day represents UTC days whether the length is 86399, 86400 or 86401 SI seconds.</para>
      /// <para>Applications should use the function iauDtf2d to convert from calendar date and time of day into 2-part quasi Julian Date, as it implements the leap-second-ambiguity convention just described.</para></description></item>
      /// <item><description>The warning status "dubious year" flags UTCs that predate the introduction of the time scale or that are too far in the future to be trusted.  See iauDat for further details.</description></item>
      /// <item><description>UT1-UTC is tabulated in IERS bulletins.  It increases by exactly one second at the end of each positive UTC leap second, introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This practice is under review, and in the future UT1-UTC may grow essentially without limit.</description></item>
      /// <item><description>The geographical coordinates are with respect to the WGS84 reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the longitude required by the present function is east-positive (i.e. right-handed), in accordance with geographical convention.</description></item>
      /// <item><description>The polar motion xp,yp can be obtained from IERS bulletins.  The values are the coordinates (in radians) of the Celestial Intermediate Pole with respect to the International Terrestrial Reference System (see IERS Conventions 2003), measured along the
      /// meridians 0 and 90 deg west respectively.  For many applications, xp and yp can be set to zero.</description></item>
      /// <item><description>If hm, the height above the ellipsoid of the observing station in meters, is not known but phpa, the pressure in hPa (=mB), is available, an adequate estimate of hm can be obtained from the expression:
      /// <p style="margin-left:25px;font-family:Lucida Conosle,Monospace"><b>hm = -29.3 * tsl * log ( phpa / 1013.25 );</b></p>
      /// <para>where tsl is the approximate sea-level air temperature in K (See Astrophysical Quantities, C.W.Allen, 3rd edition, section 52).  Similarly, if the pressure phpa is not known, it can be estimated from the height of the observing station, hm, as follows:</para>
      /// <p style="margin-left:25px;font-family:Lucida Conosle,Monospace"><b>phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );</b></p>
      /// <para>Note, however, that the refraction is nearly proportional to the pressure and that an accurate phpa value is important for precise work.</para></description></item>
      /// <item><description>The argument wl specifies the observing wavelength in micrometers.  The transition from optical to radio is assumed to occur at 100 micrometers (about 3000 GHz).</description></item>
      /// <item><description>The accuracy of the result is limited by the corrections for refraction, which use a simple A*tan(z) + B*tan^3(z) model. Providing the meteorological parameters are known accurately and
      /// there are no gross local effects, the predicted astrometric coordinates should be within 0.05 arcsec (optical) or 1 arcsec (radio) for a zenith distance of less than 70 degrees, better than 30 arcsec (optical or radio) at 85 degrees and better
      /// than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.
      ///<para>Without refraction, the complementary functions iauAtco13 and iauAtoc13 are self-consistent to better than 1 microarcsecond all over the celestial sphere.  With refraction included, consistency falls off at high zenith distances, but is still better than 0.05 arcsec at 85 degrees.</para></description></item>
      /// <item><description>It is advisable to take great care with units, as even unlikely values of the input parameters are accepted and processed in accordance with the models used.</description></item>
      /// </list>
      /// </remarks>
      public int Atoc13(string type, double ob1, double ob2, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double rc, ref double dc)
      {
         return this.ObservedToCelestial(type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref rc, ref dc);
      }

      public int ObservedToCelestial(string type, double ob1, double ob2, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double rc, ref double dc)
      {
         short RetCode = 0;

         if (Is64Bit()) {
            RetCode = Atoc1364(type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref rc, ref dc);
         }
         else {
            RetCode = Atoc1332(type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref rc, ref dc);
         }

         return Convert.ToInt32(RetCode);
      }

      /// <summary>
      ///  Observed place to CIRS using the SOFA Atoi13 function.
      /// </summary>
      /// <param name="type">type of coordinates - "R", "H" or "A" (Notes 1,2)</param>
      /// <param name="ob1">observed Az, HA or RA (radians; Az is N=0,E=90)</param>
      /// <param name="ob2">observed ZD or Dec (radians)</param>
      /// <param name="utc1">UTC as a 2-part quasi Julian Date (Notes 3,4)</param>
      /// <param name="utc2">UTC as a 2-part quasi Julian Date (Notes 3,4)</param>
      /// <param name="dut1">UT1-UTC (seconds, Note 5)</param>
      /// <param name="elong">longitude (radians, east +ve, Note 6)</param>
      /// <param name="phi">geodetic latitude (radians, Note 6)</param>
      /// <param name="hm">height above the ellipsoid (meters, Notes 6,8)</param>
      /// <param name="xp">polar motion coordinates (radians, Note 7)</param>
      /// <param name="yp">polar motion coordinates (radians, Note 7)</param>
      /// <param name="phpa">pressure at the observer (hPa = mB, Note 8)</param>
      /// <param name="tc">ambient temperature at the observer (deg C)</param>
      /// <param name="rh">relative humidity at the observer (range 0-1)</param>
      /// <param name="wl">wavelength (micrometers, Note 9)</param>
      /// <param name="ri">CIRS right ascension (CIO-based, radians)</param>
      /// <param name="di">CIRS declination (radians)</param>
      /// <returns>Status: +1 = dubious year (Note 2), 0 = OK, -1 = unacceptable date</returns>
      /// <remarks>
      /// <para>Notes:</para>
      /// <list type="number">
      /// <item><description>"Observed" Az,ZD means the position that would be seen by a perfect geodetically aligned theodolite.  (Zenith distance is used rather than altitude in order to reflect the fact that no
      /// allowance is made for depression of the horizon.)  This is related to the observed HA,Dec via the standard rotation, using the geodetic latitude (corrected for polar motion), while the
      /// observed HA and RA are related simply through the Earth rotation angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus means the position that would be seen by a perfect equatorial
      /// with its polar axis aligned to the Earth's axis of rotation.</description></item>
      /// <item><description>Only the first character of the type argument is significant. "R" or "r" indicates that ob1 and ob2 are the observed right ascension and declination;  "H" or "h" indicates that they are
      /// hour angle (west +ve) and declination;  anything else ("A" or "a" is recommended) indicates that ob1 and ob2 are azimuth (north zero, east 90 deg) and zenith distance.</description></item>
      /// <item><description>utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any convenient way between the two arguments, for example where utc1 is the Julian Day Number and utc2 is the fraction of a day.
      /// <para>However, JD cannot unambiguously represent UTC during a leap second unless special measures are taken.  The convention in the present function is that the JD day represents UTC days whether the length is 86399, 86400 or 86401 SI seconds.</para>
      /// <para>Applications should use the function iauDtf2d to convert from calendar date and time of day into 2-part quasi Julian Date, as it implements the leap-second-ambiguity convention just described.</para></description></item>
      /// <item><description>The warning status "dubious year" flags UTCs that predate the introduction of the time scale or that are too far in the future to be trusted.  See iauDat for further details.</description></item>
      /// <item><description>UT1-UTC is tabulated in IERS bulletins.  It increases by exactly one second at the end of each positive UTC leap second, introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
      /// practice is under review, and in the future UT1-UTC may grow essentially without limit.</description></item>
      /// <item><description>The geographical coordinates are with respect to the WGS84 reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the longitude required by the present function is east-positive
      /// (i.e. right-handed), in accordance with geographical convention.</description></item>
      /// <item><description>The polar motion xp,yp can be obtained from IERS bulletins.  The values are the coordinates (in radians) of the Celestial Intermediate Pole with respect to the International Terrestrial
      /// Reference System (see IERS Conventions 2003), measured along the meridians 0 and 90 deg west respectively.  For many applications, xp and yp can be set to zero.</description></item>
      /// <item><description>If hm, the height above the ellipsoid of the observing station in meters, is not known but phpa, the pressure in hPa (=mB), is available, an adequate estimate of hm can be obtained from the expression:
      /// <p style="margin-left:25px;font-family:Lucida Conosle,Monospace"><b>hm = -29.3 * tsl * log ( phpa / 1013.25 );</b></p>
      /// <para>where tsl is the approximate sea-level air temperature in K (See Astrophysical Quantities, C.W.Allen, 3rd edition, section 52).  Similarly, if the pressure phpa is not known, it can be estimated from the height of the observing station, hm, as follows:</para>
      /// <p style="margin-left:25px;font-family:Lucida Conosle,Monospace"><b>phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );</b></p>
      /// <para>Note, however, that the refraction is nearly proportional to the pressure and that an accurate phpa value is important for precise work.</para></description></item>
      /// <item><description>The argument wl specifies the observing wavelength in micrometers.  The transition from optical to radio is assumed to occur at 100 micrometers (about 3000 GHz).</description></item>
      /// <item><description>The accuracy of the result is limited by the corrections for refraction, which use a simple A*tan(z) + B*tan^3(z) model. Providing the meteorological parameters are known accurately and
      /// there are no gross local effects, the predicted astrometric coordinates should be within 0.05 arcsec (optical) or 1 arcsec (radio) for a zenith distance of less than 70 degrees, better
      /// than 30 arcsec (optical or radio) at 85 degrees and better than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.
      /// <para>Without refraction, the complementary functions iauAtio13 and iauAtoi13 are self-consistent to better than 1 microarcsecond all over the celestial sphere.  With refraction included,
      /// consistency falls off at high zenith distances, but is still better than 0.05 arcsec at 85 degrees.</para></description></item>
      /// <item><description>It is advisable to take great care with units, as even unlikely values of the input parameters are accepted and processed in accordance with the models used.</description></item>
      /// </list>
      /// </remarks>
      public int Atoi13(string type, double ob1, double ob2, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double ri, ref double di)
      {
         return this.ObservedToIntermediate(type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref ri, ref di);
      }

      public int ObservedToIntermediate(string type, double ob1, double ob2, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double ri, ref double di)
      {
         short RetCode = 0;

         if (Is64Bit()) {
            RetCode = Atoi1364(type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref ri, ref di);

         }
         else {
            RetCode = Atoi1332(type, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl, ref ri, ref di);
         }

         return Convert.ToInt32(RetCode);
      }

      /// <summary>
      /// Time scale transformation:  International Atomic Time, TAI, to Coordinated Universal Time, UTC.
      /// </summary>
      /// <param name="tai1">TAI as a 2-part Julian Date (Note 1)</param>
      /// <param name="tai2">TAI as a 2-part Julian Date (Note 1)</param>
      /// <param name="utc1">UTC as a 2-part quasi Julian Date (Notes 1-3)</param>
      /// <param name="utc2">UTC as a 2-part quasi Julian Date (Notes 1-3)</param>
      /// <returns>Status: +1 = dubious year (Note 4), 0 = OK, -1 = unacceptable date</returns>
      /// <remarks>
      /// Notes:
      /// <list type="number">
      /// <item><description>tai1+tai2 is Julian Date, apportioned in any convenient way between the two arguments, for example where tai1 is the Julian Day Number and tai2 is the fraction of a day.  The returned utc1
      /// and utc2 form an analogous pair, except that a special convention is used, to deal with the problem of leap seconds - see the next note.</description></item>
      /// <item><description>JD cannot unambiguously represent UTC during a leap second unless special measures are taken.  The convention in the present function is that the JD day represents UTC days whether the
      /// length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era there were smaller jumps (in either direction) each time the linear UTC(TAI) expression was changed, and these "mini-leaps are also included in the SOFA convention.</description></item>
      /// <item><description>The function iauD2dtf can be used to transform the UTC quasi-JD into calendar date and clock time, including UTC leap second handling.</description></item>
      /// <item><description>The warning status "dubious year" flags UTCs that predate the introduction of the time scale or that are too far in the future to be trusted.  See iauDat for further details.</description></item>
      /// </list>
      /// </remarks>
      public int TaiUtc(double tai1, double tai2, ref double utc1, ref double utc2)
      {
         short RetCode = 0;

         if (Is64Bit()) {
            RetCode = Taiutc64(tai1, tai2, ref utc1, ref utc2);
         }
         else {
            RetCode = Taiutc32(tai1, tai2, ref utc1, ref utc2);
         }

         return Convert.ToInt32(RetCode);
      }

      /// <summary>
      /// Time scale transformation:  International Atomic Time, TAI, to Terrestrial Time, TT.
      /// </summary>
      /// <param name="tai1">TAI as a 2-part Julian Date</param>
      /// <param name="tai2">TAI as a 2-part Julian Date</param>
      /// <param name="tt1">TT as a 2-part Julian Date</param>
      /// <param name="tt2">TT as a 2-part Julian Date</param>
      /// <returns>Status:  0 = OK</returns>
      /// <remarks>
      /// Notes:
      /// <list type="number">
      /// <item><description> tai1+tai2 is Julian Date, apportioned in any convenient way between the two arguments, for example where tai1 is the Julian Day Number and tai2 is the fraction of a day.  The returned
      /// tt1,tt2 follow suit.</description></item>
      /// </list>
      /// </remarks>
      public int TaiTt(double tai1, double tai2, ref double tt1, ref double tt2)
      {

         short RetCode = 0;

         if (Is64Bit()) {
            RetCode = Taitt64(tai1, tai2, ref tt1, ref tt2);
         }
         else {
            RetCode = Taitt32(tai1, tai2, ref tt1, ref tt2);
         }

         return Convert.ToInt32(RetCode);
      }

      /// <summary>
      /// Time scale transformation:  Terrestrial Time, TT, to International Atomic Time, TAI.
      /// </summary>
      /// <param name="tt1">TT as a 2-part Julian Date</param>
      /// <param name="tt2">TT as a 2-part Julian Date</param>
      /// <param name="tai1">TAI as a 2-part Julian Date</param>
      /// <param name="tai2">TAI as a 2-part Julian Date</param>
      /// <returns>Status:  0 = OK</returns>
      /// <remarks>
      /// Note
      /// <list type="number">
      /// <item><description>tt1+tt2 is Julian Date, apportioned in any convenient way between the two arguments, for example where tt1 is the Julian Day Number and tt2 is the fraction of a day.  The returned tai1,tai2 follow suit.</description></item>
      /// </list>
      /// </remarks>
      public int TtTai(double tt1, double tt2, ref double tai1, ref double tai2)
      {
         short RetCode = 0;

         if (Is64Bit()) {
            RetCode = Tttai64(tt1, tt2, ref tai1, ref tai2);
         }
         else {
            RetCode = Tttai32(tt1, tt2, ref tai1, ref tai2);
         }

         return Convert.ToInt32(RetCode);
      }

      /// <summary>
      /// Convert hours, minutes, seconds to radians.
      /// </summary>
      /// <param name="s">sign:  '-' = negative, otherwise positive</param>
      /// <param name="ihour">Hours</param>
      /// <param name="imin">Minutes</param>
      /// <param name="sec">Seconds</param>
      /// <param name="rad">Angle in radians</param>
      /// <returns>Status:  0 = OK, 1 = ihour outside range 0-23, 2 = imin outside range 0-59, 3 = sec outside range 0-59.999...</returns>
      /// <remarks>
      /// Notes:
      /// <list type="number">
      /// <item><description>The result is computed even if any of the range checks fail.</description></item>
      /// <item><description>Negative ihour, imin and/or sec produce a warning status, but the absolute value is used in the conversion.</description></item>
      /// <item><description>If there are multiple errors, the status value reflects only the first, the smallest taking precedence.</description></item>
      /// </list>
      /// </remarks>
      public int Tf2a(string s, int ihour, int imin, double sec, ref double rad)
      {

         short RetCode = 0;

         if (string.IsNullOrEmpty(s)) {
            s = " "; // Fix any invalid sign values
         }

         if (Is64Bit()) {
            Tf2a64(s.ToCharArray()[0], System.Convert.ToInt16(Convert.ToInt16(ihour)), System.Convert.ToInt16(Convert.ToInt16(imin)), sec, ref rad);
         }
         else {
            Tf2a32(s.ToCharArray()[0], System.Convert.ToInt16(Convert.ToInt16(ihour)), System.Convert.ToInt16(Convert.ToInt16(imin)), sec, ref rad);
         }

         return Convert.ToInt32(RetCode);
      }

      /// <summary>
      /// Time scale transformation:  Coordinated Universal Time, UTC, to International Atomic Time, TAI.
      /// </summary>
      /// <param name="utc1">UTC as a 2-part quasi Julian Date (Notes 1-4)</param>
      /// <param name="utc2">UTC as a 2-part quasi Julian Date (Notes 1-4)</param>
      /// <param name="tai1">TAI as a 2-part Julian Date (Note 5)</param>
      /// <param name="tai2">TAI as a 2-part Julian Date (Note 5)</param>
      /// <returns>Status: +1 = dubious year (Note 3) 0 = OK -1 = unacceptable date</returns>
      /// <remarks>
      /// Notes:
      /// <list type="number">
      /// <item><description>utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any convenient way between the two arguments, for example where utc1 is the Julian Day Number and utc2 is the fraction of a day.</description></item>
      /// <item><description>JD cannot unambiguously represent UTC during a leap second unless special measures are taken.  The convention in the present function is that the JD day represents UTC days whether the
      /// length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era there were smaller jumps (in either direction) each time the linear UTC(TAI) expression was changed, and these "mini-leaps" are also included in the SOFA convention.</description></item>
      /// <item><description>The warning status "dubious year" flags UTCs that predate the introduction of the time scale or that are too far in the future to be trusted.  See iauDat for further details.</description></item>
      /// <item><description>The function iauDtf2d converts from calendar date and time of day into 2-part Julian Date, and in the case of UTC implements the leap-second-ambiguity convention described above.</description></item>
      /// <item><description>The returned TAI1,TAI2 are such that their sum is the TAI Julian Date.</description></item>
      /// </list>
      /// </remarks>
      public int UtcTai(double utc1, double utc2, ref double tai1, ref double tai2)
      {

         short RetCode = 0;

         if (Is64Bit()) {
            Utctai64(utc1, utc2, ref tai1, ref tai2);
         }
         else {
            Utctai32(utc1, utc2, ref tai1, ref tai2);
         }

         return Convert.ToInt32(RetCode);
      }

      #endregion

      #region DLL Entry Points SOFA (32bit)
      [DllImport(SOFA32DLL, EntryPoint = "iauAf2a")] private static extern short Af2a32(char s, short ideg, short iamin, double asec, ref double rad);

      [DllImport(SOFA32DLL, EntryPoint = "iauAnp")] private static extern double Anp32(double a);

      [DllImport(SOFA32DLL, EntryPoint = "iauAtci13")] private static extern void Atci1332(double rc, double dc, double pr, double pd, double px, double rv, double date1, double date2, ref double ri, ref double di, ref double eo);

      [DllImport(SOFA32DLL, EntryPoint = "iauAtco13")] private static extern short Atco1332(double rc, double dc, double pr, double pd, double px, double rv, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double aob, ref double zob, ref double hob, ref double dob, ref double rob, ref double eo);

      [DllImport(SOFA32DLL, EntryPoint = "iauAtic13")] private static extern void Atic1332(double ri, double di, double date1, double date2, ref double rc, ref double dc, ref double eo);

      [DllImport(SOFA32DLL, EntryPoint = "iauAtoc13")] private static extern short Atoc1332(string type, double ob1, double ob2, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double rc, ref double dc);

      [DllImport(SOFA32DLL, EntryPoint = "iauAtio13")] private static extern short Atio1332(double ri, double di, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double aob, ref double zob, ref double hob, ref double dob, ref double rob);

      [DllImport(SOFA32DLL, EntryPoint = "iauAtoi13")] private static extern short Atoi1332(string type, double ob1, double ob2, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double ri, ref double di);

      [DllImport(SOFA32DLL, EntryPoint = "iauDtf2d")] private static extern short Dtf2d32(string scale, int iy, int im, int id, int ihr, int imn, double sec, ref double d1, ref double d2);

      [DllImport(SOFA32DLL, EntryPoint = "iauEo06a")] private static extern double Eo06a32(double date1, double date2);

      [DllImport(SOFA32DLL, EntryPoint = "iauTaitt")] private static extern short Taitt32(double tai1, double tai2, ref double tt1, ref double tt2);

      [DllImport(SOFA32DLL, EntryPoint = "iauTttai")] private static extern short Tttai32(double tt1, double tt2, ref double tai1, ref double tai2);

      [DllImport(SOFA32DLL, EntryPoint = "iauTf2a")] private static extern short Tf2a32(char s, short ihour, short imin, double sec, ref double rad);

      [DllImport(SOFA32DLL, EntryPoint = "iauUtctai")] private static extern short Utctai32(double utc1, double utc2, ref double tai1, ref double tai2);

      [DllImport(SOFA32DLL, EntryPoint = "iauTaiutc")] private static extern short Taiutc32(double tai1, double tai2, ref double utc1, ref double utc2);

      #endregion

      #region DLL Entry Points SOFA (64bit)
      [DllImport(SOFA64DLL, EntryPoint = "iauAf2a")] private static extern short Af2a64(char s, short ideg, short iamin, double asec, ref double rad);

      [DllImport(SOFA64DLL, EntryPoint = "iauAnp")] private static extern double Anp64(double a);

      [DllImport(SOFA64DLL, EntryPoint = "iauAtci13")] private static extern void Atci1364(double rc, double dc, double pr, double pd, double px, double rv, double date1, double date2, ref double ri, ref double di, ref double eo);

      [DllImport(SOFA64DLL, EntryPoint = "iauAtco13")] private static extern short Atco1364(double rc, double dc, double pr, double pd, double px, double rv, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double aob, ref double zob, ref double hob, ref double dob, ref double rob, ref double eo);

      [DllImport(SOFA64DLL, EntryPoint = "iauAtic13")] private static extern void Atic1364(double ri, double di, double date1, double date2, ref double rc, ref double dc, ref double eo);

      [DllImport(SOFA64DLL, EntryPoint = "iauAtoc13")] private static extern short Atoc1364(string type, double ob1, double ob2, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double rc, ref double dc);

      [DllImport(SOFA64DLL, EntryPoint = "iauAtio13")] private static extern short Atio1364(double ri, double di, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double aob, ref double zob, ref double hob, ref double dob, ref double rob);

      [DllImport(SOFA64DLL, EntryPoint = "iauAtoi13")] private static extern short Atoi1364(string type, double ob1, double ob2, double utc1, double utc2, double dut1, double elong, double phi, double hm, double xp, double yp, double phpa, double tc, double rh, double wl, ref double ri, ref double di);

      [DllImport(SOFA64DLL, EntryPoint = "iauDtf2d")] private static extern short Dtf2d64(string scale, int iy, int im, int id, int ihr, int imn, double sec, ref double d1, ref double d2);

      [DllImport(SOFA64DLL, EntryPoint = "iauEo06a")] private static extern double Eo06a64(double date1, double date2);

      [DllImport(SOFA64DLL, EntryPoint = "iauTaitt")] private static extern short Taitt64(double tai1, double tai2, ref double tt1, ref double tt2);

      [DllImport(SOFA64DLL, EntryPoint = "iauTttai")] private static extern short Tttai64(double tt1, double tt2, ref double tai1, ref double tai2);

      [DllImport(SOFA64DLL, EntryPoint = "iauTf2a")] private static extern short Tf2a64(char s, short ihour, short imin, double sec, ref double rad);

      [DllImport(SOFA64DLL, EntryPoint = "iauUtctai")] private static extern short Utctai64(double utc1, double utc2, ref double tai1, ref double tai2);

      [DllImport(SOFA64DLL, EntryPoint = "iauTaiutc")] private static extern short Taiutc64(double tai1, double tai2, ref double utc1, ref double utc2);

      #endregion

      #region Private Support Code
      //Constants for SHGetSpecialFolderPath shell call
      private const int CSIDL_PROGRAM_FILES = 38; //0x0026
      private const int CSIDL_PROGRAM_FILESX86 = 42; //0x002a,
      private const int CSIDL_WINDOWS = 36; // 0x0024,
      private const int CSIDL_PROGRAM_FILES_COMMONX86 = 44; // 0x002c,

      //DLL to provide the path to Program Files(x86)\Common Files folder location that is not avialable through the .NET framework
      /// <summary>
      /// Get path to a system folder
      /// </summary>
      /// <param name="hwndOwner">SUpply null / nothing to use "current user"</param>
      /// <param name="lpszPath">returned string folder path</param>
      /// <param name="nFolder">Folder Number from CSIDL enumeration e.g. CSIDL_PROGRAM_FILES_COMMONX86 = 44 = 0x2c</param>
      /// <param name="fCreate">Indicates whether the folder should be created if it does not already exist. If this value is nonzero,
      /// the folder is created. If this value is zero, the folder is not created</param>
      /// <returns>TRUE if successful; otherwise, FALSE.</returns>
      /// <remarks></remarks>
      [DllImport("shell32.dll")] private static extern bool SHGetSpecialFolderPath(IntPtr hwndOwner, System.Text.StringBuilder lpszPath, int nFolder, bool fCreate);

      /// <summary>
      /// Loads a library DLL
      /// </summary>
      /// <param name="lpFileName">Full path to the file to load</param>
      /// <returns>A pointer to the loaded DLL image</returns>
      /// <remarks>This is a wrapper for the Windows kernel32 function LoadLibraryA</remarks>
      [DllImport("kernel32.dll", SetLastError = true, EntryPoint = "LoadLibraryA")] public static extern IntPtr LoadLibrary(string lpFileName);

      /// <summary>
      /// Unloads a library DLL
      /// </summary>
      /// <param name="hModule">Pointer to the loaded library returned by the LoadLibrary function.</param>
      /// <returns>True or false depending on whether the library was released.</returns>
      /// <remarks></remarks>
      [DllImport("kernel32.dll", SetLastError = true, EntryPoint = "FreeLibrary")] public static extern bool FreeLibrary(IntPtr hModule);

      /// <summary>
      /// Indicates whether we are running as a 32bit or 64bit application
      /// </summary>
      /// <returns>True if the application is 64bit, False for 32bit</returns>
      /// <remarks></remarks>
      private bool Is64Bit()
      {
         if (IntPtr.Size == 8) //Check whether we are running on a 32 or 64bit system.
         {
            return true;
         }
         else {
            return false;
         }
      }

      #endregion

   }

}