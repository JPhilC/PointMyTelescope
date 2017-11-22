using ASCOM.Astrometry;
using System.Runtime.InteropServices;

namespace ASCOM.Astrometry
{
   sealed class GlobalItems
   {
      // Physical contants
      internal const double MOON_RADIUS = 1737.0; // km
      internal const double EARTH_RADIUS = 6378.0; // km
      internal const double SUN_RADIUS = 696342.0; // km
      internal const double MERCURY_RADIUS = 2439.7; // km
      internal const double VENUS_RADIUS = 2439.7; // km
      internal const double MARS_RADIUS = 3396.2; // km
      internal const double JUPITER_RADIUS = 69911.0; // km
      internal const double SATURN_RADIUS = 6051.8; // km
      internal const double NEPTUNE_RADIUS = 24767.0; // km
      internal const double URANUS_RADIUS = 24973.0; // km
      internal const double PLUTO_RADIUS = 1153.0; // km

      // Fixed event definitions
      internal const double SUN_RISE = -50.0 / 60.0; // Degrees
      internal const double CIVIL_TWIGHLIGHT = -6.0; // Degrees
      internal const double NAUTICAL_TWIGHLIGHT = -12.0; // Degrees
      internal const double AMATEUR_ASRONOMICAL_TWIGHLIGHT = -15.0; // Degrees
      internal const double ASTRONOMICAL_TWIGHLIGHT = -18.0; // Degrees

      // Conversion factors
      internal const double HOURS2DEG = 15.0;
      internal const double DEG2HOURS = 1.0 / 15.0;
      internal const double DEG2HOURSSOLSID = 1.0 / 15.04107;
      internal const double SECONDS2DAYS = 1.0 / (60.0 * 60.0 * 24.0);
      internal const double AU2KILOMETRE = 149597870.691;

      //NOVAS.COM Constants
      internal const short FN1 = 1;
      internal const short FN0 = 0;
      internal const double T0 = 2451545.0; //TDB Julian date of epoch J2000.0.
      internal const double KMAU = 149597870.0; //Astronomical Unit in kilometers.
      internal const double MAU = 149597870000.0; //Astronomical Unit in meters.
      internal const double C = 173.14463348; // Speed of light in AU/Day.
      internal const double GS = 1.3271243800000001E+20; // Heliocentric gravitational constant.
      internal const double EARTHRAD = 6378.14; //Radius of Earth in kilometers.
      internal const double F = 0.00335281; //Earth ellipsoid flattening.
      internal const double OMEGA = 0.00007292115; //Rotational angular velocity of Earth in radians/sec.
      internal const double TWOPI = 6.2831853071795862; //Value of pi in radians.
      internal const double RAD2SEC = 206264.80624709636; //Angle conversion constants.
      internal const double DEG2RAD = 0.017453292519943295;
      internal const double RAD2DEG = 57.295779513082323;

      //General constants
      internal const int TAI_UTC_OFFSET = 37; // Current TAI - UTC offset (seconds), changed to 37 at 31st December 2016
      internal const double TT_TAI_OFFSET = 32.184; //32.184 seconds
      internal const double MJDBASE = 2400000.5; //This is the offset of Modified Julian dates from true Julian dates
      internal const double SECPERDAY = 86400.0;
      internal const double J2000BASE = 2451545.0;

      internal const double RACIO_DEFAULT_VALUE = double.NaN; //NOVAS3: Default value that if still present will indicate that this value was not updated

      //Profile store names
      internal const string ASTROMETRY_SUBKEY = "Astrometry";
      internal const string UTC_TAI_OFFSET_VALUENAME = "Leap Seconds";
   }
}
   #region AstroUtil Enums and Structures
   /// <summary>
   /// Type of event for which an ephemeris is required
   /// </summary>
   /// <remarks></remarks>
   [Guid("946C9620-B292-4807-9B75-FC828AB1700B"),
      ComVisible(true)]
   public enum EventType
   {
      SunRiseSunset = 0,
      MoonRiseMoonSet = 1,
      CivilTwilight = 2,
      NauticalTwilight = 3,
      AmateurAstronomicalTwilight = 4,
      AstronomicalTwilight = 5,
      MercuryRiseSet = 6,
      VenusRiseSet = 7,
      MarsRiseSet = 8,
      JupiterRiseSet = 9,
      SaturnRiseSet = 10,
      UranusRiseSet = 11,
      NeptuneRiseSet = 12,
      PlutoRiseSet = 13
   }
   #endregion

   #region NOVAS2 Enums
   /// <summary>
   /// Type of body, Major Planet, Moon, Sun or Minor Planet
   /// </summary>
   /// <remarks></remarks>
   [Guid("A1D2C046-F7BC-474f-8D95-6E7B761DEECB"),
      ComVisible(true)]
   public enum BodyType
   {
      /// <summary>
      /// Luna
      /// </summary>
      /// <remarks></remarks>
      Moon = 0,

      /// <summary>
      /// The Sun
      /// </summary>
      /// <remarks></remarks>
      Sun = 0,

      /// <summary>
      /// Major planet
      /// </summary>
      /// <remarks></remarks>
      MajorPlanet = 0,

      /// <summary>
      /// Minor planet
      /// </summary>
      /// <remarks></remarks>
      MinorPlanet = 1,

      /// <summary>
      /// Comet
      /// </summary>
      /// <remarks></remarks>
      Comet = 2
   }

   /// <summary>
   /// Co-ordinate origin: centre of Sun or solar system barycentre
   /// </summary>
   /// <remarks></remarks>
   [Guid("9591FC6A-3EF1-41ae-9FE1-FE0C76686A85"),
      ComVisible(true)]
   public enum Origin
   {
      /// <summary>
      /// Centre of mass of the solar system
      /// </summary>
      /// <remarks></remarks>
      Barycentric = 0,
      /// <summary>
      /// Centre of mass of the Sun
      /// </summary>
      /// <remarks></remarks>
      Heliocentric = 1
   }

   /// <summary>
   /// Body number starting with Mercury = 1
   /// </summary>
   /// <remarks></remarks>
   [Guid("C839867F-E152-44e1-8356-2DB450329EDC"),
      ComVisible(true)]
   public enum Body
   {
      /// <summary>
      /// Mercury
      /// </summary>
      /// <remarks></remarks>
      Mercury = 1,
      /// <summary>
      /// Venus
      /// </summary>
      /// <remarks></remarks>
      Venus = 2,
      /// <summary>
      /// Earth
      /// </summary>
      /// <remarks></remarks>
      Earth = 3,
      /// <summary>
      /// Mars
      /// </summary>
      /// <remarks></remarks>
      Mars = 4,
      /// <summary>
      /// Jupiter
      /// </summary>
      /// <remarks></remarks>
      Jupiter = 5,
      /// <summary>
      /// Saturn
      /// </summary>
      /// <remarks></remarks>
      Saturn = 6,
      /// <summary>
      /// Uranus
      /// </summary>
      /// <remarks></remarks>
      Uranus = 7,
      /// <summary>
      /// Neptune
      /// </summary>
      /// <remarks></remarks>
      Neptune = 8,
      /// <summary>
      /// Pluto
      /// </summary>
      /// <remarks></remarks>
      Pluto = 9,
      /// <summary>
      /// Sun
      /// </summary>
      /// <remarks></remarks>
      Sun = 10,
      /// <summary>
      /// Moon
      /// </summary>
      /// <remarks></remarks>
      Moon = 11
   }

   /// <summary>
   /// Type of refraction correction
   /// </summary>
   /// <remarks></remarks>
   [Guid("32914B41-41C3-4f68-A974-E9ABCE0BA03A"),
      ComVisible(true)]
   public enum RefractionOption
   {
      /// <summary>
      /// No refraction correction will be applied
      /// </summary>
      /// <remarks></remarks>
      NoRefraction = 0,
      /// <summary>
      /// Refraction will be applied based on "standard" weather values of temperature = 10.0C and sea level pressure = 1010 millibar
      /// </summary>
      /// <remarks></remarks>
      StandardRefraction = 1,
      /// <summary>
      /// Refraction will be applied based on the temperature and pressure supplied in the site location structure
      /// </summary>
      /// <remarks></remarks>
      LocationRefraction = 2
   }

   /// <summary>
   /// Type of transformation: Epoch, Equator and Equinox or all three
   /// </summary>
   /// <remarks></remarks>
   [Guid("6ADE707E-1D7E-471a-94C9-F9FCC56755B1"),
      ComVisible(true)]
   public enum TransformationOption
   {
      /// <summary>
      /// Change epoch only
      /// </summary>
      /// <remarks></remarks>
      ChangeEpoch = 1,
      /// <summary>
      /// Change equator and equinox
      /// </summary>
      /// <remarks></remarks>
      ChangeEquatorAndEquinox = 2,
      /// <summary>
      /// Change equator, equinox and epoch
      /// </summary>
      /// <remarks></remarks>
      ChangeEquatorAndEquinoxAndEpoch = 3
   }

   /// <summary>
   /// Direction of nutation correction
   /// </summary>
   /// <remarks></remarks>
   [Guid("394E0981-3344-4cff-8ABA-E19A775AAD29"),
      ComVisible(true)]
   public enum NutationDirection
   {
      /// <summary>
      /// Convert mean equator and equinox to true equator and equinox
      /// </summary>
      /// <remarks></remarks>
      MeanToTrue = 0,
      /// <summary>
      /// Convert true equator and equinox to mean equator and equinox
      /// </summary>
      /// <remarks></remarks>
      TrueToMean = 1
   }
   #endregion

   #region NOVAS3 Enums
   /// <summary>
   /// Direction of transformation: ITRS to Terrestrial Intermediate or vice versa
   /// </summary>
   /// <remarks></remarks>
   [Guid("45EAC3DA-08FB-49E2-B852-114312933742"),
      ComVisible(true)]
   public enum TransformationDirection
   {
      ITRSToTerrestrialIntermediate = 0,
      TerrestrialIntermediateToITRS = 1
   }
   /// <summary>
   /// Location of observer
   /// </summary>
   /// <remarks></remarks>
   [Guid("8FFEAC07-F976-4fd8-8547-3DCFF25F5FA3"),
      ComVisible(true)]
   public enum ObserverLocation
   {
      /// <summary>
      /// Observer at centre of the earth
      /// </summary>
      /// <remarks></remarks>
      EarthGeoCenter = 0,
      /// <summary>
      /// Observer on earth's surface
      /// </summary>
      /// <remarks></remarks>
      EarthSurface = 1,
      /// <summary>
      /// Observer in near-earth spacecraft
      /// </summary>
      /// <remarks></remarks>
      SpaceNearEarth = 2
   }

   /// <summary>
   /// Calculation accuracy
   /// </summary>
   /// <remarks>
   /// In full-accuracy mode,
   /// <list type="bullet">
   ///<item>nutation calculations use the IAU 2000A model [iau2000a, nutation_angles];</item>
   ///<item>gravitational deflection is calculated using three bodies: Sun, Jupiter, and Saturn [grav_def];</item>
   ///<item>the equation of the equinoxes includes the entire series when computing the “complementary terms" [ee_ct];</item>
   ///<item>geocentric positions of solar system bodies are adjusted for light travel time using split, or two-part,
   /// Julian dates in calls to ephemeris and iterate with a convergence tolerance of 10-12 days [light_time, ephemeris];</item>
   ///<item>ephemeris calls the appropriate solar system ephemeris using split, or two-part, Julian dates primarily to support
   /// light-time calculations [ephemeris, solarsystem_hp, light_time].</item>
   /// </list>
   ///<para>In reduced-accuracy mode,</para>
   /// <list type="bullet">
   /// <item>nutation calculations use the 2000K model, which is the default for this mode;</item>
   /// <item>gravitational deflection is calculated using only one body, the Sun [grav_def];</item>
   /// <item>the equation of the equinoxes excludes terms smaller than 2 microarcseconds when computing the "complementary terms" [ee_ct];</item>
   /// <item>geocentric positions of solar system bodies are adjusted for light travel time using single-value Julian dates
   /// in calls to ephemeris and iterate with a convergence tolerance of 10-9 days [light-time, ephemeris, solarsystem];</item>
   /// <item>ephemeris calls the appropriate solar system ephemeris using single-value Julian dates [ephemeris, solarsystem].</item>
   /// </list>
   /// <para>In full-accuracy mode, the IAU 2000A nutation series (1,365 terms) is used [iau2000a]. Evaluating the series for nutation is
   /// usually the main computational burden in NOVAS, so using reduced-accuracy mode improves execution time, often noticeably.
   /// In reduced-accuracy mode, the NOVAS 2000K nutation series (488 terms) is used by default [nu2000k]. This mode can be used
   /// when the accuracy requirements are not better than 0.1 milliarcsecond for stars or 3.5 milliarcseconds for solar system bodies.
   /// Selecting this approach can reduce the time required for Earth-rotation computations by about two-thirds.</para>
   /// </remarks>
   [Guid("F10B748F-4F90-4acf-9EB0-76D50293E9A9"),
      ComVisible(true)]
   public enum Accuracy
   {
      /// <summary>
      /// Full accuracy
      /// </summary>
      /// <remarks>Suitable when precision of better than 0.1 milliarcsecond for stars or 3.5 milliarcseconds for solar system bodies is required.</remarks>
      Full = 0, //... full accuracy
                /// <summary>
                /// Reduced accuracy
                /// </summary>
                /// <remarks>Suitable when precision of less than 0.1 milliarcsecond for stars or 3.5 milliarcseconds for solar system bodies is required.</remarks>
      Reduced = 1 //... reduced accuracy
   }

   /// <summary>
   /// Coordinate system of the output position
   /// </summary>
   /// <remarks>Used by function Place</remarks>
   [Guid("0EF9BC38-B790-4416-8FEF-E03758B6B630"),
      ComVisible(true)]
   public enum CoordSys
   {
      /// <summary>
      /// GCRS or "local GCRS"
      /// </summary>
      /// <remarks></remarks>
      GCRS = 0,
      /// <summary>
      /// True equator and equinox of date
      /// </summary>
      /// <remarks></remarks>
      EquinoxOfDate = 1,
      /// <summary>
      /// True equator and CIO of date
      /// </summary>
      /// <remarks></remarks>
      CIOOfDate = 2,
      /// <summary>
      /// Astrometric coordinates, i.e., without light deflection or aberration.
      /// </summary>
      /// <remarks></remarks>
      Astrometric = 3
   }

   /// <summary>
   /// Type of sidereal time
   /// </summary>
   /// <remarks></remarks>
   [Guid("7722AE51-F475-4c69-8B35-B2EDBD297C66"),
      ComVisible(true)]
   public enum GstType
   {
      /// <summary>
      /// Greenwich mean sidereal time
      /// </summary>
      /// <remarks></remarks>
      GreenwichMeanSiderealTime = 0,
      /// <summary>
      /// Greenwich apparent sidereal time
      /// </summary>
      /// <remarks></remarks>
      GreenwichApparentSiderealTime = 1
   }

   /// <summary>
   /// Computation method
   /// </summary>
   /// <remarks></remarks>
   [Guid("8D9E6EF5-CE9C-4ba9-8B24-C0FA5067D8FA"),
      ComVisible(true)]
   public enum Method
   {
      /// <summary>
      /// Based on CIO
      /// </summary>
      /// <remarks></remarks>
      CIOBased = 0,
      /// <summary>
      /// Based on equinox
      /// </summary>
      /// <remarks></remarks>
      EquinoxBased = 1
   }

   /// <summary>
   /// Output vector reference system
   /// </summary>
   /// <remarks></remarks>
   [Guid("CD7AEAC0-1BFA-447e-A43E-62C231B0FC55"),
      ComVisible(true)]
   public enum OutputVectorOption
   {
      /// <summary>
      /// Referred to GCRS axes
      /// </summary>
      /// <remarks></remarks>
      ReferredToGCRSAxes = 0,
      /// <summary>
      /// Referred to the equator and equinox of date
      /// </summary>
      /// <remarks></remarks>
      ReferredToEquatorAndEquinoxOfDate = 1
   }

   /// <summary>
   /// Type of pole ofset
   /// </summary>
   /// <remarks>Used by CelPole.</remarks>
   [Guid("AF69D7CC-A59C-4fcc-BE17-C2F568957BFD"),
      ComVisible(true)]
   public enum PoleOffsetCorrection
   {
      /// <summary>
      /// For corrections to angular coordinates of modeled pole referred to mean ecliptic of date, that is, delta-delta-psi
      /// and delta-delta-epsilon.
      /// </summary>
      /// <remarks></remarks>
      ReferredToMeanEclipticOfDate = 1,
      /// <summary>
      /// For corrections to components of modeled pole unit vector referred to GCRS axes, that is, dx and dy.
      /// </summary>
      /// <remarks></remarks>
      ReferredToGCRSAxes = 2
   }

   /// <summary>
   /// Direction of frame conversion
   /// </summary>
   /// <remarks>Used by FrameTie method.</remarks>
   [Guid("3AC3E32A-EDCE-4234-AA50-CDB346851C5D"),
      ComVisible(true)]
   public enum FrameConversionDirection
   {
      /// <summary>
      /// Dynamical to ICRS transformation.
      /// </summary>
      /// <remarks></remarks>
      DynamicalToICRS = -1,
      /// <summary>
      /// ICRS to dynamical transformation.
      /// </summary>
      /// <remarks></remarks>
      ICRSToDynamical = 1
   }

   /// <summary>
   /// Location of observer, determining whether the gravitational deflection due to the earth itself is applied.
   /// </summary>
   /// <remarks>Used by GravDef method.</remarks>
   [Guid("C39A798C-53F7-460b-853F-DA5389B4324D"),
      ComVisible(true)]
   public enum EarthDeflection
   {
      /// <summary>
      /// No earth deflection (normally means observer is at geocenter)
      /// </summary>
      /// <remarks></remarks>
      NoEarthDeflection = 0,
      /// <summary>
      /// Add in earth deflection (normally means observer is on or above surface of earth, including earth orbit)
      /// </summary>
      /// <remarks></remarks>
      AddEarthDeflection = 1
   }

   /// <summary>
   /// Reference system in which right ascension is given
   /// </summary>
   /// <remarks></remarks>
   [Guid("DF215CCF-2C25-48e5-A357-A8300C1EA027"),
      ComVisible(true)]
   public enum ReferenceSystem
   {
      /// <summary>
      /// GCRS
      /// </summary>
      /// <remarks></remarks>
      GCRS = 1,
      /// <summary>
      /// True equator and equinox of date
      /// </summary>
      /// <remarks></remarks>
      TrueEquatorAndEquinoxOfDate = 2
   }

   /// <summary>
   /// Type of equinox
   /// </summary>
   /// <remarks></remarks>
   [Guid("5EDEE8B3-E223-4fb7-924F-4C56E8373380"),
      ComVisible(true)]
   public enum EquinoxType
   {
      /// <summary>
      /// Mean equinox
      /// </summary>
      /// <remarks></remarks>
      MeanEquinox = 0,
      /// <summary>
      /// True equinox
      /// </summary>
      /// <remarks></remarks>
      TrueEquinox = 1
   }

   /// <summary>
   /// Type of transformation
   /// </summary>
   /// <remarks></remarks>
   [Guid("8BBA934E-D874-48a2-A3E2-C842A7FFFB35"),
      ComVisible(true)]
   public enum TransformationOption3
   {
      /// <summary>
      /// Change epoch only
      /// </summary>
      /// <remarks></remarks>
      ChangeEpoch = 1,
      /// <summary>
      /// Change equator and equinox; sane epoch
      /// </summary>
      /// <remarks></remarks>
      ChangeEquatorAndEquinox = 2,
      /// <summary>
      /// Change equator, equinox and epoch
      /// </summary>
      /// <remarks></remarks>
      ChangeEquatorAndEquinoxAndEpoch = 3,
      /// <summary>
      /// change equator and equinox J2000.0 to ICRS
      /// </summary>
      /// <remarks></remarks>
      ChangeEquatorAndEquinoxJ2000ToICRS = 4,
      /// <summary>
      /// change ICRS to equator and equinox of J2000.0
      /// </summary>
      /// <remarks></remarks>
      ChangeICRSToEquatorAndEquinoxOfJ2000 = 5
   }

   /// <summary>
   /// Type of object
   /// </summary>
   /// <remarks></remarks>
   [Guid("5BBB931B-358C-40ac-921C-B48373F01348"),
      ComVisible(true)]
   public enum ObjectType
   {
      /// <summary>
      /// Major planet, sun or moon
      /// </summary>
      /// <remarks></remarks>
      MajorPlanetSunOrMoon = 0,
      /// <summary>
      /// Minor planet
      /// </summary>
      /// <remarks></remarks>
      MinorPlanet = 1,
      /// <summary>
      /// Object located outside the solar system
      /// </summary>
      /// <remarks></remarks>
      ObjectLocatedOutsideSolarSystem = 2
   }

   /// <summary>
   /// Body or location
   /// </summary>
   /// <remarks>This numbering convention is used by ephemeris routines; do not confuse with the Body enum, which is used in most
   /// other places within NOVAS3.
   /// <para>
   /// The numbering convention for 'target' and'center' is:
   /// <pre>
   ///             0  =  Mercury           7 = Neptune
   ///             1  =  Venus             8 = Pluto
   ///             2  =  Earth             9 = Moon
   ///             3  =  Mars             10 = Sun
   ///             4  =  Jupiter          11 = Solar system bary.
   ///             5  =  Saturn           12 = Earth-Moon bary.
   ///             6  =  Uranus           13 = Nutations (long. and obliq.)</pre>
   /// </para>
   /// <para>
   /// If nutations are desired, set 'target' = 14; 'center' will be ignored on that call.
   /// </para>
   ///</remarks>
   [Guid("60E342F9-3CC3-4b98-8045-D61B5A7D974B"),
      ComVisible(true)]
   public enum Target
   {
      /// <summary>
      /// Mercury
      /// </summary>
      /// <remarks></remarks>
      Mercury = 0,
      /// <summary>
      /// Venus
      /// </summary>
      /// <remarks></remarks>
      Venus = 1,
      /// <summary>
      /// Earth
      /// </summary>
      /// <remarks></remarks>
      Earth = 2,
      /// <summary>
      /// Mars
      /// </summary>
      /// <remarks></remarks>
      Mars = 3,
      /// <summary>
      /// Jupiter
      /// </summary>
      /// <remarks></remarks>
      Jupiter = 4,
      /// <summary>
      /// Saturn
      /// </summary>
      /// <remarks></remarks>
      Saturn = 5,
      /// <summary>
      /// Uranus
      /// </summary>
      /// <remarks></remarks>
      Uranus = 6,
      /// <summary>
      /// Neptune
      /// </summary>
      /// <remarks></remarks>
      Neptune = 7,
      /// <summary>
      /// Pluto
      /// </summary>
      /// <remarks></remarks>
      Pluto = 8,
      /// <summary>
      /// Moon
      /// </summary>
      /// <remarks></remarks>
      Moon = 9,
      /// <summary>
      /// Sun
      /// </summary>
      /// <remarks></remarks>
      Sun = 10,
      /// <summary>
      /// Solar system barycentre
      /// </summary>
      /// <remarks></remarks>
      SolarSystemBarycentre = 11,
      /// <summary>
      /// Earth moon barycentre
      /// </summary>
      /// <remarks></remarks>
      EarthMoonBarycentre = 12,
      /// <summary>
      /// Nutations
      /// </summary>
      /// <remarks></remarks>
      Nutations = 13
   }
   #endregion

   #region Public NOVAS 2 Structures
   /// <summary>
   /// Structure to hold body type, number and name
   /// </summary>
   /// <remarks>Designates a celestial object.
   /// </remarks>
   [Guid("558F644F-E112-4e88-9D79-20063BB25C3E"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential)]
   public struct BodyDescription
   {
      /// <summary>
      /// Type of body
      /// </summary>
      /// <remarks>
      /// 0 = Major planet, Sun, or Moon
      /// 1 = Minor planet
      /// </remarks>
      public BodyType Type;
      /// <summary>
      /// body number
      /// </summary>
      /// <remarks><pre>
      /// For 'type' = 0: Mercury = 1, ..., Pluto = 9, Sun = 10, Moon = 11
      /// For 'type' = 1: minor planet number
      /// </pre></remarks>
      public Body Number;
      /// <summary>
      /// Name of the body (limited to 99 characters)
      /// </summary>
      /// <remarks></remarks>
      [MarshalAsAttribute(UnmanagedType.BStr, SizeConst = 100)] public string Name; //char[100]
   }

   /// <summary>
   /// Structure to hold astrometric catalogue data
   /// </summary>
   /// <remarks>
   /// The astrometric catalog data for a star; equator and equinox and units will depend on the catalog.
   /// While this structure can be used as a generic container for catalog data, all high-level
   /// NOVAS-C functions require J2000.0 catalog data with FK5-type units (shown in square brackets below).
   /// </remarks>
   [Guid("6320FEDA-8582-4048-988A-7D4DE7978C71"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential, CharSet = CharSet.Ansi)]
   public struct CatEntry
   {
      /// <summary>
      /// 3-character catalog designator.
      /// </summary>
      /// <remarks></remarks>
      [MarshalAsAttribute(UnmanagedType.BStr, SizeConst = 4)] public string Catalog; //char[4] was <MarshalAsAttribute(UnmanagedType.ByValTStr, SizeConst:=4)> Was this before changing for COM compatibility

      /// <summary>
      /// Name of star.
      /// </summary>
      /// <remarks></remarks>
      [MarshalAsAttribute(UnmanagedType.BStr, SizeConst = 51)] public string StarName; //char[51] was <MarshalAsAttribute(UnmanagedType.ByValTStr, SizeConst:=51)>
                                                                                       /// <summary>
                                                                                       /// Integer identifier assigned to star.
                                                                                       /// </summary>
                                                                                       /// <remarks></remarks>
      public int StarNumber;

      /// <summary>
      /// Mean right ascension [hours].
      /// </summary>
      /// <remarks></remarks>
      public double RA;

      /// <summary>
      /// Mean declination [degrees].
      /// </summary>
      /// <remarks></remarks>
      public double Dec;

      /// <summary>
      /// Proper motion in RA [seconds of time per century].
      /// </summary>
      /// <remarks></remarks>
      public double ProMoRA;

      /// <summary>
      /// Proper motion in declination [arcseconds per century].
      /// </summary>
      /// <remarks></remarks>
      public double ProMoDec;

      /// <summary>
      /// Parallax [arcseconds].
      /// </summary>
      /// <remarks></remarks>
      public double Parallax;

      /// <summary>
      /// Radial velocity [kilometers per second]
      /// </summary>
      /// <remarks></remarks>
      public double RadialVelocity;
   }

   /// <summary>
   /// Structure to hold site information
   /// </summary>
   /// <remarks>
   /// Data for the observer's location.  The atmospheric parameters are used only by the refraction
   /// function called from function 'equ_to_hor'. Additional parameters can be added to this
   /// structure if a more sophisticated refraction model is employed.
   /// </remarks>
   [Guid("ED02B64A-320F-47cd-90D9-3DF2DF07602D"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential, CharSet = CharSet.Ansi)]
   public struct SiteInfo
   {
      /// <summary>
      /// Geodetic latitude in degrees; north positive.
      /// </summary>
      /// <remarks></remarks>
      public double Latitude; //geodetic latitude in degrees; north positive.
                              /// <summary>
                              /// Geodetic longitude in degrees; east positive.
                              /// </summary>
                              /// <remarks></remarks>
      public double Longitude; //geodetic longitude in degrees; east positive.
                               /// <summary>
                               /// Height of the observer in meters.
                               /// </summary>
                               /// <remarks></remarks>
      public double Height; //height of the observer in meters.
                            /// <summary>
                            /// Temperature (degrees Celsius).
                            /// </summary>
                            /// <remarks></remarks>
      public double Temperature; //temperature (degrees Celsius).
                                 /// <summary>
                                 /// Atmospheric pressure (millibars)
                                 /// </summary>
                                 /// <remarks></remarks>
      public double Pressure; //atmospheric pressure (millibars)
   }

   /// <summary>
   /// Structure to hold a position vector
   /// </summary>
   /// <remarks>Object position vector
   /// </remarks>
   [Guid("69651C90-75F5-4f46-8D0F-22D186151D45"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential, CharSet = CharSet.Ansi)]
   public struct PosVector
   {
      /// <summary>
      /// x co-ordinate
      /// </summary>
      /// <remarks></remarks>
      public double x;
      /// <summary>
      /// y co-ordinate
      /// </summary>
      /// <remarks></remarks>
      public double y;
      /// <summary>
      /// z co-ordinate
      /// </summary>
      /// <remarks></remarks>
      public double z;
   }

   /// <summary>
   /// Structure to hold a velocity vector
   /// </summary>
   /// <remarks>Object velocity vector
   /// </remarks>
   [Guid("F18240B0-00CC-4ff7-9A94-AC835387F959"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential, CharSet = CharSet.Ansi)]
   public struct VelVector
   {
      /// <summary>
      /// x velocity component
      /// </summary>
      /// <remarks></remarks>
      public double x;
      /// <summary>
      /// y velocity component
      /// </summary>
      /// <remarks></remarks>
      public double y;
      /// <summary>
      /// z velocity component
      /// </summary>
      /// <remarks></remarks>
      public double z;
   }

   /// <summary>
   /// Structure to hold Sun and Moon fundamental arguments
   /// </summary>
   /// <remarks>Fundamental arguments, in radians
   ///</remarks>
   [Guid("5EE28FFB-39CD-4d23-BF62-11EE4C581681"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential, CharSet = CharSet.Ansi)]
   public struct FundamentalArgs
   {
      /// <summary>
      /// l (mean anomaly of the Moon)
      /// </summary>
      /// <remarks></remarks>
      public double l;
      /// <summary>
      /// l' (mean anomaly of the Sun)
      /// </summary>
      /// <remarks></remarks>
      public double ldash;
      /// <summary>
      /// F (L - omega; L = mean longitude of the Moon)
      /// </summary>
      /// <remarks></remarks>
      public double F;
      /// <summary>
      /// D (mean elongation of the Moon from the Sun)
      /// </summary>
      /// <remarks></remarks>
      public double D;
      /// <summary>
      /// Omega (mean longitude of the Moon's ascending node)
      /// </summary>
      /// <remarks></remarks>
      public double Omega;
   }
   #endregion

   #region Public NOVAS 3 Structures

   /// <summary>
   /// Catalogue entry structure
   /// </summary>
   /// <remarks>Basic astrometric data for any celestial object located outside the solar system; the catalog data for a star.
   /// <para>This structure is identical to the NOVAS2 CatEntry structure expect that, for some reason, the StarName and Catalog fields
   /// have been swapped in the NOVAS3 structure.</para>
   /// <para>
   /// Please note that some units have changed from those used in NOVAS2 as follows:
   /// <list type="bullet">
   /// <item>proper motion in right ascension: from seconds per century to milliarcseconds per year</item>
   /// <item>proper motion in declination: from arcseconds per century to milliarcseconds per year</item>
   /// <item>parallax: from arcseconds to milliarcseconds</item>
   /// </list>
   /// </para>
   /// </remarks>
   [Guid("5325E96C-BD24-4470-A0F6-E917B05805E1"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential, CharSet = CharSet.Ansi)]
   public struct CatEntry3
   {
      /// <summary>
      /// Name of celestial object. (maximum 50 characters)
      /// </summary>
      /// <remarks></remarks>
      [MarshalAsAttribute(UnmanagedType.ByValTStr, SizeConst = 51)] public string StarName;

      /// <summary>
      /// 3-character catalog designator.
      /// </summary>
      /// <remarks></remarks>
      [MarshalAsAttribute(UnmanagedType.ByValTStr, SizeConst = 4)] public string Catalog;

      /// <summary>
      /// Integer identifier assigned to object.
      /// </summary>
      /// <remarks></remarks>
      public int StarNumber;

      /// <summary>
      /// ICRS right ascension (hours)
      /// </summary>
      /// <remarks></remarks>
      public double RA;

      /// <summary>
      /// ICRS declination (degrees)
      /// </summary>
      /// <remarks></remarks>
      public double Dec;

      /// <summary>
      /// ICRS proper motion in right ascension (milliarcseconds/year)
      /// </summary>
      /// <remarks></remarks>
      public double ProMoRA;

      /// <summary>
      /// ICRS proper motion in declination (milliarcseconds/year)
      /// </summary>
      /// <remarks></remarks>
      public double ProMoDec;

      /// <summary>
      /// Parallax (milli-arcseconds)
      /// </summary>
      /// <remarks></remarks>
      public double Parallax;

      /// <summary>
      /// Radial velocity (km/s)
      /// </summary>
      /// <remarks></remarks>
      public double RadialVelocity;
   }

   /// <summary>
   /// Celestial object structure
   /// </summary>
   /// <remarks>Designates a celestial object</remarks>
   [Guid("AEFE0EA0-D013-46a9-B77D-6D0FDD661005"),
      ComVisible(true)]
   public struct Object3
   {
      /// <summary>
      /// Type of object
      /// </summary>
      /// <remarks></remarks>
      public ObjectType Type;
      /// <summary>
      /// Object identification number
      /// </summary>
      /// <remarks></remarks>
      public Body Number;
      /// <summary>
      /// Name of object(maximum 50 characters)
      /// </summary>
      /// <remarks></remarks>
      public string Name;
      /// <summary>
      /// Catalogue entry for the object
      /// </summary>
      /// <remarks></remarks>
      public CatEntry3 Star;
   }

   /// <summary>
   /// Celestial object's place in the sky
   /// </summary>
   /// <remarks></remarks>
   [Guid("9AD852C3-A895-4f69-AEC0-C9CA44283FA0"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential, CharSet = CharSet.Ansi)]
   public struct SkyPos
   {
      /// <summary>
      /// Unit vector toward object (dimensionless)
      /// </summary>
      /// <remarks></remarks>
      [MarshalAsAttribute(UnmanagedType.ByValArray, SizeConst = 3, ArraySubType = UnmanagedType.R8)] public double[] RHat;
      /// <summary>
      /// Apparent, topocentric, or astrometric right ascension (hours)
      /// </summary>
      /// <remarks></remarks>
      public double RA;
      /// <summary>
      /// Apparent, topocentric, or astrometric declination (degrees)
      /// </summary>
      /// <remarks></remarks>
      public double Dec;
      /// <summary>
      /// True (geometric, Euclidian) distance to solar system body or 0.0 for star (AU)
      /// </summary>
      /// <remarks></remarks>
      public double Dis;
      /// <summary>
      /// Radial velocity (km/s)
      /// </summary>
      /// <remarks></remarks>
      public double RV;
   }

   /// <summary>
   /// Observer’s position and velocity in a near-Earth spacecraft.
   /// </summary>
   /// <remarks></remarks>
   [Guid("15737EA5-E4FA-40da-8BDA-B8CF96D89E43"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential)]
   public struct InSpace
   {
      /// <summary>
      /// Geocentric position vector (x, y, z), components in km with respect to true equator and equinox of date
      /// </summary>
      /// <remarks></remarks>
      [MarshalAsAttribute(UnmanagedType.ByValArray, SizeConst = 3, ArraySubType = UnmanagedType.R8)] public double[] ScPos;
      /// <summary>
      /// Geocentric velocity vector (x_dot, y_dot, z_dot), components in km/s with respect to true equator and equinox of date
      /// </summary>
      /// <remarks></remarks>
      [MarshalAsAttribute(UnmanagedType.ByValArray, SizeConst = 3, ArraySubType = UnmanagedType.R8)] public double[] ScVel;
   }

   /// <summary>
   /// Right ascension of the Celestial Intermediate Origin (CIO) with respect to the GCRS.
   /// </summary>
   /// <remarks></remarks>
   [Guid("4959930F-0CDB-4324-A0E0-F60A351454B7"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential)]
   public struct RAOfCio
   {
      /// <summary>
      /// TDB Julian date
      /// </summary>
      /// <remarks></remarks>
      public double JdTdb;
      /// <summary>
      /// Right ascension of the CIO with respect to the GCRS (arcseconds)
      /// </summary>
      /// <remarks></remarks>
      public double RACio;
   }

   /// <summary>
   /// Parameters of observer's location
   /// </summary>
   /// <remarks>This structure is identical to the NOVAS2 SiteInfo structure but is included so that NOVAS3 naming
   /// conventions are maintained, making it easier to relate this code to the NOVAS3 documentation and C code.</remarks>
   [Guid("277380A5-6599-448f-9232-1C280073D3CD"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential)]
   public struct OnSurface
   {
      /// <summary>
      /// Geodetic (ITRS) latitude; north positive (degrees)
      /// </summary>
      /// <remarks></remarks>
      public double Latitude;
      /// <summary>
      /// Geodetic (ITRS) longitude; east positive (degrees)
      /// </summary>
      /// <remarks></remarks>
      public double Longitude;
      /// <summary>
      /// Observer's height above sea level
      /// </summary>
      /// <remarks></remarks>
      public double Height;
      /// <summary>
      /// Observer's location's ambient temperature (degrees Celsius)
      /// </summary>
      /// <remarks></remarks>
      public double Temperature;
      /// <summary>
      /// Observer's location's atmospheric pressure (millibars)
      /// </summary>
      /// <remarks></remarks>
      public double Pressure;
   }

   /// <summary>
   /// General specification for the observer's location
   /// </summary>
   /// <remarks></remarks>
   [Guid("64A25FDD-3687-45e0-BEAF-18C361E5E340"),
      ComVisible(true),
      StructLayoutAttribute(LayoutKind.Sequential)]
   public struct Observer
   {
      /// <summary>
      /// Code specifying the location of the observer: 0=at geocenter; 1=surface of earth; 2=near-earth spacecraft
      /// </summary>
      /// <remarks></remarks>
      public ObserverLocation Where;
      /// <summary>
      /// Data for an observer's location on the surface of the Earth (where = 1)
      /// </summary>
      /// <remarks></remarks>
      public OnSurface OnSurf;
      /// <summary>
      /// Data for an observer's location on a near-Earth spacecraft (where = 2)
      /// </summary>
      /// <remarks></remarks>
      public InSpace NearEarth;
   }
   #endregion

   #region Internal NOVAS3 Structures
   [StructLayoutAttribute(LayoutKind.Sequential, CharSet = CharSet.Ansi)]
   internal struct JDHighPrecision
   {
      public double JDPart1;
      public double JDPart2;
   }

   //Internal version of Object3 with correct marshalling hints and type for Number field
   [StructLayoutAttribute(LayoutKind.Sequential, CharSet = CharSet.Ansi)]
   internal struct Object3Internal
   {
      public ObjectType Type;
      public short Number;
      [MarshalAsAttribute(UnmanagedType.ByValTStr, SizeConst = 51)] public string Name;
      public CatEntry3 Star;
   }

   [StructLayoutAttribute(LayoutKind.Sequential)]
   internal struct RAOfCioArray
   {

      internal RAOfCio Value1;
      internal RAOfCio Value2;
      internal RAOfCio Value3;
      internal RAOfCio Value4;
      internal RAOfCio Value5;
      internal RAOfCio Value6;
      internal RAOfCio Value7;
      internal RAOfCio Value8;
      internal RAOfCio Value9;
      internal RAOfCio Value10;
      internal RAOfCio Value11;
      internal RAOfCio Value12;
      internal RAOfCio Value13;
      internal RAOfCio Value14;
      internal RAOfCio Value15;
      internal RAOfCio Value16;
      internal RAOfCio Value17;
      internal RAOfCio Value18;
      internal RAOfCio Value19;
      internal RAOfCio Value20;

      internal void Initialise()
      {
         Value1.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value2.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value3.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value4.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value5.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value6.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value7.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value8.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value9.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value10.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value11.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value12.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value13.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value14.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value15.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value16.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value17.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value18.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value19.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
         Value20.RACio = GlobalItems.RACIO_DEFAULT_VALUE;
      }
   }
   #endregion
