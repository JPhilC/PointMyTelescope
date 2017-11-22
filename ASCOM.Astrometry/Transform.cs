using System;
using System.Globalization;
using System.Runtime.InteropServices;

//Transform component implementation


namespace ASCOM.Astrometry.Transform
{
   /// <summary>
   /// Coordinate transform component; J2000 - apparent - local topocentric
   /// </summary>
   /// <remarks>Use this component to transform between J2000, apparent and local topocentric (JNow) coordinates or
   /// vice versa. To use the component, instantiate it, then use one of SetJ2000 or SetJNow or SetApparent to
   /// initialise with known values. Now use the RAJ2000, DECJ200, RAJNow, DECJNow, RAApparent and DECApparent etc.
   /// properties to read off the required transformed values.
   ///<para>The component can be reused simply by setting new co-ordinates with a Set command, there
   /// is no need to create a new component each time a transform is required.</para>
   /// <para>Transforms are effected through the ASCOM NOVAS.Net engine that encapsulates the USNO NOVAS 3.1 library.
   /// The USNO NOVAS reference web page is:
   /// <href>http://www.usno.navy.mil/USNO/astronomical-applications/software-products/novas</href>
   /// and the NOVAS 3.1 user guide is included in the ASCOM Developer Components install.
   /// </para>
   /// </remarks>
   public class Transform:IDisposable 
   {
      private bool disposedValue = false; // To detect redundant calls
      SOFA.SOFA SOFA;
      private double RAJ2000Value;
      private double RATopoValue;
      private double DECJ2000Value;
      private double DECTopoValue;
      private double SiteElevValue;
      private double SiteLatValue;
      private double SiteLongValue;
      private double SiteTempValue;
      private double RAApparentValue;
      private double DECApparentValue;
      private double AzimuthTopoValue;
      private double ElevationTopoValue;
      private double JulianDateTTValue;
      private double JulianDateUTCValue;
      private bool RefracValue;
      private bool RequiresRecalculate;
      private SetBy LastSetBy;

      private const double HOURS2RADIANS = Math.PI / 12.0;
      private const double DEGREES2RADIANS = Math.PI / 180.0;
      private const double RADIANS2HOURS = 12.0 / Math.PI;
      private const double RADIANS2DEGREES = 180.0 / Math.PI;
      private const double TWOPI = 2.0 * Math.PI;

      private const string DATE_FORMAT = "dd/MM/yyyy HH:mm:ss.fff";

      private enum SetBy
      {
         Never,
         J2000,
         Apparent,
         Topocentric,
         AzimuthElevation,
         Refresh
      }

      #region New and IDisposable
      public Transform()
      {

         SOFA = new SOFA.SOFA();

         RAJ2000Value = System.Convert.ToDouble(double.NaN); //Initialise to invalid values in case these are read before they are set
         DECJ2000Value = System.Convert.ToDouble(double.NaN);
         RATopoValue = System.Convert.ToDouble(double.NaN);
         DECTopoValue = System.Convert.ToDouble(double.NaN);
         SiteElevValue = System.Convert.ToDouble(double.NaN);
         SiteLatValue = System.Convert.ToDouble(double.NaN);
         SiteLongValue = System.Convert.ToDouble(double.NaN);
         RefracValue = false;
         LastSetBy = SetBy.Never;
         RequiresRecalculate = true;
         JulianDateTTValue = 0; // Initialise to a value that forces the current PC date time to be used in determining the TT Julian date of interest
         CheckGAC();
      }


      // IDisposable
      protected virtual void Dispose(bool disposing)
      {
         if (!this.disposedValue) {
            if (disposing) {
               // Free other state (managed objects).

               if (SOFA != null) {
                  SOFA.Dispose();
                  SOFA = null;
               }
            }
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


      #region ITransform Implementation
      /// <summary>
      /// Gets or sets the site latitude
      /// </summary>
      /// <value>Site latitude</value>
      /// <returns>Latitude in degrees</returns>
      /// <remarks>Positive numbers north of the equator, negative numbers south.</remarks>
      public double SiteLatitude
      {
         get
         {
            CheckSet("SiteLatitude", SiteLatValue, "Site latitude has not been set");
            return SiteLatValue;
         }
         set
         {
            if (SiteLatValue != value) {
               RequiresRecalculate = true;
            }
            SiteLatValue = value;
         }
      }

      /// <summary>
      /// Gets or sets the site longitude
      /// </summary>
      /// <value>Site longitude</value>
      /// <returns>Longitude in degrees</returns>
      /// <remarks>Positive numbers east of the Greenwich meridian, negative numbes west of the Greenwich meridian.</remarks>
      public double SiteLongitude
      {
         get
         {
            CheckSet("SiteLongitude", SiteLongValue, "Site longitude has not been set");
            return SiteLongValue;
         }
         set
         {
            if (SiteLongValue != value) {
               RequiresRecalculate = true;
            }
            SiteLongValue = value;
         }
      }

      /// <summary>
      /// Gets or sets the site elevation above sea level
      /// </summary>
      /// <value>Site elevation</value>
      /// <returns>Elevation in metres</returns>
      /// <remarks></remarks>
      public double SiteElevation
      {
         get
         {
            CheckSet("SiteElevation", SiteElevValue, "Site elevation has not been set");
            return SiteElevValue;
         }
         set
         {
            if (SiteElevValue != value) {
               RequiresRecalculate = true;
            }
            SiteElevValue = value;
         }
      }

      /// <summary>
      /// Gets or sets the site ambient temperature
      /// </summary>
      /// <value>Site ambient temperature</value>
      /// <returns>Temperature in degrees Celsius</returns>
      /// <remarks></remarks>
      public double SiteTemperature
      {
         get
         {
            CheckSet("SiteTemperature", SiteTempValue, "Site temperature has not been set");
            return SiteTempValue;
         }
         set
         {
            if (SiteTempValue != value) {
               RequiresRecalculate = true;
            }
            SiteTempValue = value;
         }
      }

      /// <summary>
      /// Gets or sets a flag indicating whether refraction is calculated for topocentric co-ordinates
      /// </summary>
      /// <value>True / false flag indicating refaction is included / omitted from topocentric co-ordinates</value>
      /// <returns>Boolean flag</returns>
      /// <remarks></remarks>
      public bool Refraction
      {
         get
         {
            return RefracValue;
         }
         set
         {
            if (RefracValue != value) {
               RequiresRecalculate = true;
            }
            RefracValue = value;
         }
      }

      /// <summary>
      /// Causes the transform component to recalculate values derrived from the last Set command
      /// </summary>
      /// <remarks>Use this when you have set J2000 co-ordinates and wish to ensure that the mount points to the same
      /// co-ordinates allowing for local effects that change with time such as refraction.
      /// <para><b style="color:red">Note:</b> As of Platform 6 SP2 use of this method is not required, refresh is always performed automatically when required.</para></remarks>
      public void Refresh()
      {
         Recalculate();
      }

      /// <summary>
      /// Sets the known J2000 Right Ascension and Declination coordinates that are to be transformed
      /// </summary>
      /// <param name="RA">RA in J2000 co-ordinates</param>
      /// <param name="DEC">DEC in J2000 co-ordinates</param>
      /// <remarks></remarks>
      public void SetJ2000(double RA, double DEC)
      {
         LastSetBy = SetBy.J2000;
         if ((RA != RAJ2000Value) || (DEC != DECJ2000Value)) {
            RequiresRecalculate = true;
         }
         RAJ2000Value = ValidateRA("SetJ2000", RA);
         DECJ2000Value = ValidateDec("SetJ2000", DEC);
      }

      /// <summary>
      /// Sets the known apparent Right Ascension and Declination coordinates that are to be transformed
      /// </summary>
      /// <param name="RA">RA in apparent co-ordinates</param>
      /// <param name="DEC">DEC in apparent co-ordinates</param>
      /// <remarks></remarks>
      public void SetApparent(double RA, double DEC)
      {
         LastSetBy = SetBy.Apparent;
         if ((RA != RAApparentValue) || (DEC != DECApparentValue)) {
            RequiresRecalculate = true;
         }
         RAApparentValue = ValidateRA("SetApparent", RA);
         DECApparentValue = ValidateDec("SetApparent", DEC);
      }

      ///<summary>
      /// Sets the known local topocentric Right Ascension and Declination coordinates that are to be transformed
      /// </summary>
      /// <param name="RA">RA in local topocentric co-ordinates</param>
      /// <param name="DEC">DEC in local topocentric co-ordinates</param>
      /// <remarks></remarks>
      public void SetTopocentric(double RA, double DEC)
      {
         LastSetBy = SetBy.Topocentric;
         if ((RA != RATopoValue) || (DEC != DECTopoValue)) {
            RequiresRecalculate = true;
         }
         RATopoValue = ValidateRA("SetTopocentric", RA);
         DECTopoValue = ValidateDec("SetTopocentric", DEC);
      }

      /// <summary>
      /// Sets the topocentric azimuth and elevation
      /// </summary>
      /// <param name="Azimuth">Topocentric Azimuth in degrees</param>
      /// <param name="Elevation">Topocentric elevation in degrees</param>
      /// <remarks></remarks>
      public void SetAzimuthElevation(double Azimuth, double Elevation)
      {
         LastSetBy = SetBy.AzimuthElevation;
         RequiresRecalculate = true;
         AzimuthTopoValue = Azimuth;
         ElevationTopoValue = Elevation;
      }

      /// <summary>
      /// Returns the Right Ascension in J2000 co-ordinates
      /// </summary>
      /// <value>J2000 Right Ascension</value>
      /// <returns>Right Ascension in hours</returns>
      /// <exception cref="Exceptions.TransformUninitialisedException">Exception thrown if an attempt is made
      /// to read a value before any of the Set methods has been used or if the value can not be derived from the
      /// information in the last Set method used. E.g. topocentric values will be unavailable if the last Set was
      /// a SetApparent and one of the Site properties has not been set.</exception>
      /// <remarks></remarks>
      ///
      public double RAJ2000
      {
         get
         {
            if (LastSetBy == SetBy.Never) {
               throw (new Exceptions.TransformUninitialisedException("Attempt to read RAJ2000 before a SetXX method has been called"));
            }
            Recalculate();
            CheckSet("RAJ2000", RAJ2000Value, "RA J2000 can not be derived from the information provided. Are site parameters set?");
            return RAJ2000Value;
         }
      }

      /// <summary>
      /// Returns the Declination in J2000 co-ordinates
      /// </summary>
      /// <value>J2000 Declination</value>
      /// <returns>Declination in degrees</returns>
      /// <exception cref="Exceptions.TransformUninitialisedException">Exception thrown if an attempt is made
      /// to read a value before any of the Set methods has been used or if the value can not be derived from the
      /// information in the last Set method used. E.g. topocentric values will be unavailable if the last Set was
      /// a SetApparent and one of the Site properties has not been set.</exception>
      /// <remarks></remarks>
      public double DECJ2000
      {
         get
         {
            return this.DecJ2000;
         }
      }

      public double DecJ2000
      {
         get
         {
            if (LastSetBy == SetBy.Never) {
               throw (new Exceptions.TransformUninitialisedException("Attempt to read DECJ2000 before a SetXX method has been called"));
            }
            Recalculate();
            CheckSet("DecJ2000", DECJ2000Value, "DEC J2000 can not be derived from the information provided. Are site parameters set?");
            return DECJ2000Value;
         }
      }

      /// <summary>
      /// Returns the Right Ascension in local topocentric co-ordinates
      /// </summary>
      /// <value>Local topocentric Right Ascension</value>
      /// <returns>Right Ascension in hours</returns>
      /// <exception cref="Exceptions.TransformUninitialisedException">Exception thrown if an attempt is made
      /// to read a value before any of the Set methods has been used or if the value can not be derived from the
      /// information in the last Set method used. E.g. topocentric values will be unavailable if the last Set was
      /// a SetApparent and one of the Site properties has not been set.</exception>
      /// <remarks></remarks>
      public double RATopocentric
      {
         get
         {
            if (LastSetBy == SetBy.Never) {
               throw (new Exceptions.TransformUninitialisedException("Attempt to read RATopocentric before a SetXX method  has been called"));
            }
            Recalculate();
            CheckSet("RATopocentric", RATopoValue, "RA topocentric can not be derived from the information provided. Are site parameters set?");
            return RATopoValue;
         }
      }

      /// <summary>
      /// Returns the Declination in local topocentric co-ordinates
      /// </summary>
      /// <value>Local topocentric Declination</value>
      /// <returns>Declination in degrees</returns>
      /// <exception cref="Exceptions.TransformUninitialisedException">Exception thrown if an attempt is made
      /// to read a value before any of the Set methods has been used or if the value can not be derived from the
      /// information in the last Set method used. E.g. topocentric values will be unavailable if the last Set was
      /// a SetApparent and one of the Site properties has not been set.</exception>
      /// <remarks></remarks>
      public double DECTopocentric
      {
         get
         {
            if (LastSetBy == SetBy.Never) {
               throw (new Exceptions.TransformUninitialisedException("Attempt to read DECTopocentric before a SetXX method has been called"));
            }
            Recalculate();
            CheckSet("DECTopocentric", DECTopoValue, "DEC topocentric can not be derived from the information provided. Are site parameters set?");
            return DECTopoValue;
         }
      }

      /// <summary>
      /// Returns the Right Ascension in apparent co-ordinates
      /// </summary>
      /// <value>Apparent Right Ascension</value>
      /// <returns>Right Ascension in hours</returns>
      /// <exception cref="Exceptions.TransformUninitialisedException">Exception thrown if an attempt is made
      /// to read a value before any of the Set methods has been used or if the value can not be derived from the
      /// information in the last Set method used. E.g. topocentric values will be unavailable if the last Set was
      /// a SetApparent and one of the Site properties has not been set.</exception>
      /// <remarks></remarks>
      public double RAApparent
      {
         get
         {
            if (LastSetBy == SetBy.Never) {
               throw (new Exceptions.TransformUninitialisedException("Attempt to read DECApparent before a SetXX method has been called"));
            }
            Recalculate();
            return RAApparentValue;
         }
      }

      /// <summary>
      /// Returns the Declination in apparent co-ordinates
      /// </summary>
      /// <value>Apparent Declination</value>
      /// <returns>Declination in degrees</returns>
      /// <exception cref="Exceptions.TransformUninitialisedException">Exception thrown if an attempt is made
      /// to read a value before any of the Set methods has been used or if the value can not be derived from the
      /// information in the last Set method used. E.g. topocentric values will be unavailable if the last Set was
      /// a SetApparent and one of the Site properties has not been set.</exception>
      /// <remarks></remarks>
      public double DECApparent
      {
         get
         {
            if (LastSetBy == SetBy.Never) {
               throw (new Exceptions.TransformUninitialisedException("Attempt to read DECApparent before a SetXX method has been called"));
            }
            Recalculate();
            return DECApparentValue;
         }
      }

      /// <summary>
      /// Returns the topocentric azimth angle of the target
      /// </summary>
      /// <value>Topocentric azimuth angle</value>
      /// <returns>Azimuth angle in degrees</returns>
      /// <exception cref="Exceptions.TransformUninitialisedException">Exception thrown if an attempt is made
      /// to read a value before any of the Set methods has been used or if the value can not be derived from the
      /// information in the last Set method used. E.g. topocentric values will be unavailable if the last Set was
      /// a SetApparent and one of the Site properties has not been set.</exception>
      /// <remarks></remarks>
      public double AzimuthTopocentric
      {
         get
         {
            if (LastSetBy == SetBy.Never) {
               throw (new Exceptions.TransformUninitialisedException("Attempt to read AzimuthTopocentric before a SetXX method has been called"));
            }
            RequiresRecalculate = true; //Force a recalculation of Azimuth
            Recalculate();
            CheckSet("AzimuthTopocentric", AzimuthTopoValue, "Azimuth topocentric can not be derived from the information provided. Are site parameters set?");
            return AzimuthTopoValue;
         }
      }

      /// <summary>
      /// Returns the topocentric elevation of the target
      /// </summary>
      /// <value>Topocentric elevation angle</value>
      /// <returns>Elevation angle in degrees</returns>
      /// <exception cref="Exceptions.TransformUninitialisedException">Exception thrown if an attempt is made
      /// to read a value before any of the Set methods has been used or if the value can not be derived from the
      /// information in the last Set method used. E.g. topocentric values will be unavailable if the last Set was
      /// a SetApparent and one of the Site properties has not been set.</exception>
      /// <remarks></remarks>
      public double ElevationTopocentric
      {
         get
         {
            if (LastSetBy == SetBy.Never) {
               throw (new Exceptions.TransformUninitialisedException("Attempt to read ElevationTopocentric before a SetXX method has been called"));
            }
            RequiresRecalculate = true; //Force a recalculation of Elevation
            Recalculate();
            CheckSet("ElevationTopocentric", ElevationTopoValue, "Elevation topocentric can not be derived from the information provided. Are site parameters set?");
            return ElevationTopoValue;
         }
      }

      /// <summary>
      /// Sets or returns the Julian date on the Terrestrial Time timescale for which the transform will be made
      /// </summary>
      /// <value>Julian date (Terrestrial Time) of the transform</value>
      /// <returns>Terrestrial Time Julian date that will be used by Transform or zero if the PC's current clock value will be used to calculate the Julian date.</returns>
      /// <remarks>This method was introduced in May 2012. Previously, Transform used the current date-time of the PC when calculating transforms;
      /// this remains the default behaviour for backward compatibility.
      /// The inital value of this parameter is 0.0, which is a special value that forces Transform to replicate original behaviour by determining the
      /// Julian date from the PC's current date and time. If this property is non zero, that particular terrestrial time Julian date is used in preference
      /// to the value derrived from the PC's clock.
      /// <para>Only one of JulianDateTT or JulianDateUTC needs to be set. Use whichever is more readily available, there is no
      /// need to set both values. Transform will use the last set value of either JulianDateTT or JulianDateUTC as the basis for its calculations.</para></remarks>
      public double JulianDateTT
      {
         get
         {
            return JulianDateTTValue;
         }
         set
         {
            double tai1 = 0;
            double tai2 = 0;
            double utc1 = 0;
            double utc2 = 0;
            JulianDateTTValue = value;
            RequiresRecalculate = true; // Force a recalculation because the Julian date has changed

            if (JulianDateTTValue != 0.0) {
               //Calculate UTC
               JulianDateUTCValue = utc1 + utc2;
            }
            else // Handle special case of 0.0
            {
               JulianDateUTCValue = 0.0;
            }
         }
      }

      /// <summary>
      /// Sets or returns the Julian date on the UTC timescale for which the transform will be made
      /// </summary>
      /// <value>Julian date (UTC) of the transform</value>
      /// <returns>UTC Julian date that will be used by Transform or zero if the PC's current clock value will be used to calculate the Julian date.</returns>
      /// <remarks>Introduced in April 2014 as an alternative to JulianDateTT. Only one of JulianDateTT or JulianDateUTC needs to be set. Use whichever is more readily available, there is no
      /// need to set both values. Transform will use the last set value of either JulianDateTT or JulianDateUTC as the basis for its calculations.</remarks>
      public double JulianDateUTC
      {
         get
         {
            return JulianDateUTCValue;
         }
         set
         {
            double tai1 = 0;
            double tai2 = 0;
            double tt1 = 0;
            double tt2 = 0;
            JulianDateUTCValue = value;
            RequiresRecalculate = true; // Force a recalculation because the Julian date has changed

            if (JulianDateUTCValue != 0.0) {
               // Calculate Terrestrial Time equivalent
               JulianDateTTValue = tt1 + tt2;
            }
            else // Handle special case of 0.0
            {
               JulianDateTTValue = 0.0;
            }
         }
      }
      #endregion

      #region Support Code

      private void CheckSet(string Caller, double Value, string ErrMsg)
      {
         if (double.IsNaN(Value)) {
            throw (new Exceptions.TransformUninitialisedException(ErrMsg));
         }
      }

      private void J2000ToTopo()
      {
         double DeltaT;
         double DUT1 = 0;
         double JDUTCSofa = 0;
         double aob = 0;
         double zob = 0;
         double hob = 0;
         double dob = 0;
         double rob = 0;
         double eo = 0;
         DateTime JDUTCSofaDateTime;

         if (double.IsNaN(SiteElevValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site elevation has not been set"));
         }
         if (double.IsNaN(SiteLatValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site latitude has not been set"));
         }
         if (double.IsNaN(SiteLongValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site longitude has not been set"));
         }
         if (double.IsNaN(SiteTempValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site temperature has not been set"));
         }


         JDUTCSofa = GetJDUTCSofa();
         DeltaT = System.Convert.ToDouble(DeltatCode.DeltaTCalc(JDUTCSofa));
         DUT1 = System.Convert.ToDouble(DeltaUT(JDUTCSofa));

         JDUTCSofaDateTime = Julian2DateTime(JDUTCSofa);


         if (RefracValue) // Include refraction
         {
            SOFA.CelestialToObserved(RAJ2000Value * HOURS2RADIANS, DECJ2000Value * DEGREES2RADIANS, 0.0, 0.0, 0.0, 0.0, JDUTCSofa, 0.0, DUT1, SiteLongValue * DEGREES2RADIANS, SiteLatValue * DEGREES2RADIANS, SiteElevValue, 0.0, 0.0, 1000.0, SiteTempValue, 0.8, 0.57, ref aob, ref zob, ref hob, ref dob, ref rob, ref eo);
         }
         else // No refraction
         {
            SOFA.CelestialToObserved(RAJ2000Value * HOURS2RADIANS, DECJ2000Value * DEGREES2RADIANS, 0.0, 0.0, 0.0, 0.0, JDUTCSofa, 0.0, DUT1, SiteLongValue * DEGREES2RADIANS, SiteLatValue * DEGREES2RADIANS, SiteElevValue, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ref aob, ref zob, ref hob, ref dob, ref rob, ref eo);
         }

         RATopoValue = System.Convert.ToDouble(SOFA.Anp(rob - eo) * RADIANS2HOURS); // // Convert CIO RA to equinox of date RA by subtracting the equation of the origins and convert from radians to hours
         DECTopoValue = dob * RADIANS2DEGREES; // Convert Dec from radians to degrees
         AzimuthTopoValue = aob * RADIANS2DEGREES;
         ElevationTopoValue = System.Convert.ToDouble(90.0 - zob * RADIANS2DEGREES);

      }

      private void J2000ToApparent()
      {
         double ri = 0;
         double di = 0;
         double eo = 0;
         double JDTTSofa = 0;

         JDTTSofa = GetJDTTSofa();

         SOFA.CelestialToIntermediate(RAJ2000Value * HOURS2RADIANS, DECJ2000Value * DEGREES2RADIANS, 0.0, 0.0, 0.0, 0.0, JDTTSofa, 0.0, ref ri, ref di, ref eo);
         RAApparentValue = System.Convert.ToDouble(SOFA.Anp(ri - eo) * RADIANS2HOURS); // // Convert CIO RA to equinox of date RA by subtracting the equation of the origins and convert from radians to hours
         DECApparentValue = di * RADIANS2DEGREES; // Convert Dec from radians to degrees


      }

      private void TopoToJ2000()
      {
         double RACelestrial = 0;
         double DecCelestial = 0;
         double JDTTSofa = 0;
         double JDUTCSofa = 0;
         double DUT1 = 0;
         int RetCode;
         double aob = 0;
         double zob = 0;
         double hob = 0;
         double dob = 0;
         double rob = 0;
         double eo = 0;

         if (double.IsNaN(SiteElevValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site elevation has not been set"));
         }
         if (double.IsNaN(SiteLatValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site latitude has not been set"));
         }
         if (double.IsNaN(SiteLongValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site longitude has not been set"));
         }
         if (double.IsNaN(SiteTempValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site temperature has not been set"));
         }


         JDUTCSofa = GetJDUTCSofa();
         JDTTSofa = GetJDTTSofa();
         DUT1 = System.Convert.ToDouble(DeltaUT(JDUTCSofa));

         if (RefracValue) // Refraction is requuired
         {
            RetCode = System.Convert.ToInt32(SOFA.ObservedToCelestial("R", SOFA.Anp(RATopoValue * HOURS2RADIANS + SOFA.Eo06a(JDTTSofa, 0.0)), DECTopoValue * DEGREES2RADIANS, JDUTCSofa, 0.0, DUT1, SiteLongValue * DEGREES2RADIANS, SiteLatValue * DEGREES2RADIANS, SiteElevValue, 0.0, 0.0, 1000, SiteTempValue, 0.85, 0.57, ref RACelestrial, ref DecCelestial));
         }
         else {
            RetCode = System.Convert.ToInt32(SOFA.ObservedToCelestial("R", SOFA.Anp(RATopoValue * HOURS2RADIANS + SOFA.Eo06a(JDTTSofa, 0.0)), DECTopoValue * DEGREES2RADIANS, JDUTCSofa, 0.0, DUT1, SiteLongValue * DEGREES2RADIANS, SiteLatValue * DEGREES2RADIANS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ref RACelestrial, ref DecCelestial));
         }

         RAJ2000Value = RACelestrial * RADIANS2HOURS;
         DECJ2000Value = DecCelestial * RADIANS2DEGREES;

         // Now calculate the corresponding AzEl values from the J2000 values
         if (RefracValue) // Include refraction
         {
            SOFA.CelestialToObserved(RAJ2000Value * HOURS2RADIANS, DECJ2000Value * DEGREES2RADIANS, 0.0, 0.0, 0.0, 0.0, JDUTCSofa, 0.0, DUT1, SiteLongValue * DEGREES2RADIANS, SiteLatValue *
               DEGREES2RADIANS, SiteElevValue, 0.0, 0.0, 1000.0, SiteTempValue, 0.8, 0.57, ref aob, ref zob, ref hob, ref dob, ref rob, ref eo);
         }
         else // No refraction
         {
            SOFA.CelestialToObserved(RAJ2000Value * HOURS2RADIANS, DECJ2000Value * DEGREES2RADIANS, 0.0, 0.0, 0.0, 0.0, JDUTCSofa, 0.0, DUT1, SiteLongValue * DEGREES2RADIANS, SiteLatValue * DEGREES2RADIANS, SiteElevValue, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ref aob, ref zob, ref hob, ref dob, ref rob, ref eo);
         }

         AzimuthTopoValue = aob * RADIANS2DEGREES;
         ElevationTopoValue = System.Convert.ToDouble(90.0 - zob * RADIANS2DEGREES);


      }

      private void ApparentToJ2000()
      {
         double JulianDateTTSofa = 0;
         double RACelestial = 0;
         double DecCelestial = 0;
         double JulianDateUTCSofa = 0;
         double eo = 0;


         JulianDateTTSofa = GetJDTTSofa();
         JulianDateUTCSofa = GetJDUTCSofa();

         SOFA.IntermediateToCelestial(SOFA.Anp(RAApparentValue * HOURS2RADIANS + SOFA.Eo06a(JulianDateUTCSofa, 0.0)), DECApparentValue * DEGREES2RADIANS, JulianDateTTSofa, 0.0, ref RACelestial, ref DecCelestial, ref eo);
         RAJ2000Value = RACelestial * RADIANS2HOURS;
         DECJ2000Value = DecCelestial * RADIANS2DEGREES;

      }

      private void Recalculate() //Calculate values for derrived co-ordinates
      {
         if (RequiresRecalculate || (RefracValue == true)) {
            switch (LastSetBy) {
               case SetBy.J2000: //J2000 coordinates have bee set so calculate apparent and topocentric coords
                  //Check whether required topo values have been set
                  if ((!double.IsNaN(SiteLatValue)) && (!double.IsNaN(SiteLongValue)) && (!double.IsNaN(SiteElevValue)) && (!double.IsNaN(SiteTempValue))) {
                     J2000ToTopo(); //All required site values present so calc Topo values
                  }
                  else //Set to NaN
                  {
                     RATopoValue = System.Convert.ToDouble(double.NaN);
                     DECTopoValue = System.Convert.ToDouble(double.NaN);
                     AzimuthTopoValue = System.Convert.ToDouble(double.NaN);
                     ElevationTopoValue = System.Convert.ToDouble(double.NaN);
                  }
                  J2000ToApparent();
                  break;
               case SetBy.Topocentric: //Topocentric co-ordinates have been set so calculate J2000 and apparent coords
                  //Check whether required topo values have been set
                  if ((!double.IsNaN(SiteLatValue)) && (!double.IsNaN(SiteLongValue)) && (!double.IsNaN(SiteElevValue)) && (!double.IsNaN(SiteTempValue))) //They have so calculate remaining values
                  {
                     TopoToJ2000();
                     J2000ToApparent();
                  }
                  else //Set the topo and apaprent values to NaN
                  {
                     RAJ2000Value = System.Convert.ToDouble(double.NaN);
                     DECJ2000Value = System.Convert.ToDouble(double.NaN);
                     RAApparentValue = System.Convert.ToDouble(double.NaN);
                     DECApparentValue = System.Convert.ToDouble(double.NaN);
                     AzimuthTopoValue = System.Convert.ToDouble(double.NaN);
                     ElevationTopoValue = System.Convert.ToDouble(double.NaN);
                  }
                  break;
               case SetBy.Apparent: //Apparent values have been set so calculate J2000 values and topo values if appropriate
                  ApparentToJ2000(); //Calculate J2000 value
                                     //Check whether required topo values have been set
                  if ((!double.IsNaN(SiteLatValue)) && (!double.IsNaN(SiteLongValue)) && (!double.IsNaN(SiteElevValue)) && (!double.IsNaN(SiteTempValue))) {
                     J2000ToTopo(); //All required site values present so calc Topo values
                  }
                  else {
                     RATopoValue = System.Convert.ToDouble(double.NaN);
                     DECTopoValue = System.Convert.ToDouble(double.NaN);
                     AzimuthTopoValue = System.Convert.ToDouble(double.NaN);
                     ElevationTopoValue = System.Convert.ToDouble(double.NaN);
                  }
                  break;
               case SetBy.AzimuthElevation:
                  if ((!double.IsNaN(SiteLatValue)) && (!double.IsNaN(SiteLongValue)) && (!double.IsNaN(SiteElevValue)) && (!double.IsNaN(SiteTempValue))) {
                     AzElToJ2000();
                     J2000ToTopo();
                     J2000ToApparent();
                  }
                  else {
                     RAJ2000Value = System.Convert.ToDouble(double.NaN);
                     DECJ2000Value = System.Convert.ToDouble(double.NaN);
                     RAApparentValue = System.Convert.ToDouble(double.NaN);
                     DECApparentValue = System.Convert.ToDouble(double.NaN);
                     RATopoValue = System.Convert.ToDouble(double.NaN);
                     DECTopoValue = System.Convert.ToDouble(double.NaN);
                  }
                  break;
               default: //Neither SetJ2000 nor SetTopocentric nor SetApparent have been called, so throw an exception
                  throw (new Exceptions.TransformUninitialisedException("Can\'t recalculate Transform object values because neither SetJ2000 nor SetTopocentric nor SetApparent have been called"));
                  break;
            }
            RequiresRecalculate = false; //Reset the recalculate flag
         }
         else {
         }
      }

      private void AzElToJ2000()
      {
         int RetCode;
         double JulianDateUTCSofa = 0;
         double JulianDateTTSofa;
         double RACelestial = 0;
         double DecCelestial = 0;
         double DUT1 = 0;


         if (double.IsNaN(SiteElevValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site elevation has not been set"));
         }
         if (double.IsNaN(SiteLatValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site latitude has not been set"));
         }
         if (double.IsNaN(SiteLongValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site longitude has not been set"));
         }
         if (double.IsNaN(SiteTempValue)) {
            throw (new Exceptions.TransformUninitialisedException("Site temperature has not been set"));
         }

         JulianDateUTCSofa = GetJDUTCSofa();
         JulianDateTTSofa = GetJDTTSofa();
         DUT1 = System.Convert.ToDouble(DeltaUT(JulianDateUTCSofa));

         if (RefracValue) // Refraction is requuired
         {
            RetCode = System.Convert.ToInt32(SOFA.ObservedToCelestial("A", AzimuthTopoValue * DEGREES2RADIANS, (90.0 - ElevationTopoValue) * DEGREES2RADIANS, JulianDateUTCSofa, 0.0, DUT1, SiteLongValue * DEGREES2RADIANS, SiteLatValue * DEGREES2RADIANS, SiteElevValue, 0.0, 0.0, 1000, SiteTempValue, 0.85, 0.57, ref RACelestial, ref DecCelestial));
         }
         else {
            RetCode = System.Convert.ToInt32(SOFA.ObservedToCelestial("A", AzimuthTopoValue * DEGREES2RADIANS, (90.0 - ElevationTopoValue) * DEGREES2RADIANS, JulianDateUTCSofa, 0.0, DUT1, SiteLongValue * DEGREES2RADIANS, SiteLatValue * DEGREES2RADIANS, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ref RACelestial, ref DecCelestial));
         }

         RAJ2000Value = RACelestial * RADIANS2HOURS;
         DECJ2000Value = DecCelestial * RADIANS2DEGREES;
      }

      private double GetJDUTCSofa()
      {
         double Retval = 0;
         double utc1 = 0;
         double utc2 = 0;
         DateTime Now = default(DateTime);

         if (JulianDateUTCValue == 0.0) {
            Now = DateTime.UtcNow;
            Retval = utc1 + utc2;
         }
         else {
            Retval = JulianDateUTCValue;
         }
         return Retval;
      }

      private double GetJDTTSofa()
      {
         double Retval = 0;
         double utc1 = 0;
         double utc2 = 0;
         double tai1 = 0;
         double tai2 = 0;
         double tt1 = 0;
         double tt2 = 0;
         DateTime Now = default(DateTime);

         if (JulianDateTTValue == 0.0) {
            Now = DateTime.UtcNow;
            Retval = tt1 + tt2;
         }
         else {
            Retval = JulianDateTTValue;
         }
         return Retval;
      }

      private void CheckGAC()
      {
         string strPath = "";
         strPath = System.Convert.ToString(System.IO.Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location));
      }

      private double ValidateRA(string Caller, double RA)
      {
         if ((RA < 0.0) || (RA >= 24.0)) {
            throw (new ArgumentOutOfRangeException("RA", RA, "RA must lie between 0 and 23.9999"));
         }
         return RA;
      }

      private double ValidateDec(string Caller, double Dec)
      {
         if ((Dec < -90.0) || (Dec > 90.0)) {
            throw (new ArgumentOutOfRangeException("Dec", Dec, "Dec must lie between -90.0 to 90.0"));
         }
         return Dec;
      }

      private DateTime Julian2DateTime(double m_JulianDate)
      {
         long L = 0;
         long N = 0;
         long I = 0;
         long J = 0;
         long JDLong = 0;
         double JDFraction = 0;
         double Remainder = 0;
         int Day = 0;
         int Month = 0;
         int Year = 0;
         int Hours = 0;
         int Minutes = 0;
         int Seconds = 0;
         int MilliSeconds = 0;
         DateTime Retval = default(DateTime);

         try {
            if (m_JulianDate > 2378507.5) // 1/1/1800
            {
               JDLong = System.Convert.ToInt64(Math.Floor(m_JulianDate));
               JDFraction = m_JulianDate - Math.Floor(m_JulianDate);

               L = System.Convert.ToInt64(JDLong + 68569);
               N = System.Convert.ToInt64((4 * L) / 146097);
               L = System.Convert.ToInt64(L - System.Convert.ToInt32((System.Convert.ToInt32(146097 * N) + 3) / 4));
               I = System.Convert.ToInt64(4000 * System.Convert.ToInt32(L + 1) / 1461001);
               L = System.Convert.ToInt64(System.Convert.ToInt32(L - (System.Convert.ToInt32((1461 * I) / 4))) + 31);
               J = System.Convert.ToInt64((80 * L) / 2447);
               Day = System.Convert.ToInt32(L - System.Convert.ToInt32((2447 * J) / 80));
               L = System.Convert.ToInt64(J / 11);
               Month = System.Convert.ToInt32(System.Convert.ToDouble(J + 2) - System.Convert.ToDouble(12 * L));
               Year = System.Convert.ToInt32(System.Convert.ToInt32(100 * System.Convert.ToInt32(N - 49)) + I + L);


               JDFraction += 5.0 / (24.0 * 60.0 * 60.0 * 10000.0);

               // Allow for Julian days to start at 12:00 rather than 00:00
               if (JDFraction >= 0.5) // After midnight so add 1 to the julian day and remove half a day from the day fraction
               {
                  Day++;
                  JDFraction -= 0.5;
               }
               else // Before midnight so just add half a day
               {
                  JDFraction += 0.5;
               }

               Hours = System.Convert.ToInt32(JDFraction * 24.0);
               Remainder = System.Convert.ToDouble((JDFraction * 24.0) - Hours);

               Minutes = System.Convert.ToInt32(Remainder * 60.0);
               Remainder = System.Convert.ToDouble((Remainder * 60.0) - Minutes);

               Seconds = System.Convert.ToInt32(Remainder * 60.0);
               Remainder = System.Convert.ToDouble((Remainder * 60.0) - Seconds);

               MilliSeconds = System.Convert.ToInt32(Remainder * 1000.0);

               Retval = new DateTime(Year, Month, Day, Hours, Minutes, Seconds, MilliSeconds);
            }
            else // Early or invalid julian date so return a default value
                     {
               Retval = new DateTime(1800, 1, 10); // Return this as a default bad value
            }
         }


         catch (Exception ex) {
            Retval = new DateTime(1900, 1, 10); // Return this as a default bad value
         }

         return (Retval);
      }

      public double DeltaUT(double JulianDate)
      {
         double DUT1 = 0;

         DUT1 = System.Convert.ToDouble(GlobalItems.TAI_UTC_OFFSET + GlobalItems.TT_TAI_OFFSET - DeltatCode.DeltaTCalc(JulianDate));
         return DUT1;
      }

      #endregion

   }
}