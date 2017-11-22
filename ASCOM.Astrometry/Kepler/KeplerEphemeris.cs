using ASCOM.Utilities;
using ASCOM.Astrometry.NOVASCOM;
using System.Runtime.InteropServices;

//Kepler component implementation


namespace ASCOM.Astrometry.Kepler
{

   /// <summary>
   /// KEPLER: Ephemeris Object
   /// </summary>
   /// <remarks>
   /// The Kepler Ephemeris object contains an orbit engine which takes the orbital parameters of a solar system
   /// body, plus a a terrestrial date/time, and produces the heliocentric equatorial position and
   /// velocity vectors of the body in Cartesian coordinates. Orbital parameters are not required for
   /// the major planets, Kepler contains an ephemeris generator for these bodies that is within 0.05
   /// arc seconds of the JPL DE404 over a wide range of times, Perturbations from major planets are applied
   /// to ephemerides for minor planets.
   /// <para>The results are passed back as an array containing the two vectors.
   /// Note that this is the format expected for the ephemeris generator used by the NOVAS-COM vector
   /// astrometry engine. For more information see the description of Ephemeris.GetPositionAndVelocity().</para>
   /// <para>
   /// <b>Ephemeris Calculations</b><br />
   /// The ephemeris calculations in Kepler draw heavily from the work of
   /// Stephen Moshier moshier@world.std.com. kepler is released as a free software package, further
   /// extending the work of Mr. Moshier.</para>
   /// <para>Kepler does not integrate orbits to the current epoch. If you want the accuracy resulting from
   /// an integrated orbit, you must integrate separately and supply Kepler with elements of the current
   /// epoch. Orbit integration is on the list of things for the next major version.</para>
   /// <para>Kepler uses polynomial approximations for the major planet ephemerides. The tables
   /// of coefficients were derived by a least squares fit of periodic terms to JPL's DE404 ephemerides.
   /// The periodic frequencies used were determined by spectral analysis and comparison with VSOP87 and
   /// other analytical planetary theories. The least squares fit to DE404 covers the interval from -3000
   /// to +3000 for the outer planets, and -1350 to +3000 for the inner planets. For details on the
   /// accuracy of the major planet ephemerides, see the Accuracy Tables page. </para>
   /// <para>
   /// <b>Date and Time Systems</b><br /><br />
   /// For a detailed explanation of astronomical timekeeping systems, see A Time Tutorial on the NASA
   /// Goddard Spaceflight Center site, and the USNO Systems of Time site.
   /// <br /><br /><i>ActiveX Date values </i><br />
   /// These are the Windows standard "date serial" numbers, and are expressed in local time or
   /// UTC (see below). The fractional part of these numbers represents time within a day.
   /// They are used throughout applications such as Excel, Visual Basic, VBScript, and other
   /// ActiveX capable environments.
   /// <br /><br /><i>Julian dates </i><br />
   /// These are standard Julian "date serial" numbers, and are expressed in UTC time or Terrestrial
   /// time. The fractional part of these numbers represents time within a day. The standard ActiveX
   /// "Double" precision of 15 digits gives a resolution of about one millisecond in a full Julian date.
   /// This is sufficient for the purposes of this program.
   /// <br /><br /><i>Hourly Time Values </i><br />
   /// These are typically used to represent sidereal time and right ascension. They are simple real
   /// numbers in units of hours.
   /// <br /><br /><i>UTC Time Scale </i><br />
   /// Most of the ASCOM methods and properties that accept date/time values (either Date or Julian)
   /// assume that the date/time is in Coordinated Universal Time (UTC). Where necessary, this time
   /// is converted internally to other scales. Note that UTC seconds are based on the Cesium atom,
   /// not planetary motions. In order to keep UTC in sync with planetary motion, leap seconds are
   /// inserted periodically. The error is at most 900 milliseconds.
   /// <br /><br /><i>UT1 Time Scale </i><br />
   /// The UT1 time scale is the planetary equivalent of UTC. It it runs smoothly and varies a bit
   /// with time, but it is never more than 900 milliseconds different from UTC.
   /// <br /><br /><i>TT Time Scale </i><br />
   /// The Terrestrial Dynamical Time (TT) scale is used in solar system orbital calculations.
   /// It is based completely on planetary motions; you can think of the solar system as a giant
   /// TT clock. It differs from UT1 by an amount called "delta-t", which slowly increases with time,
   /// and is about 60 seconds right now (2001). </para>
   /// </remarks>
   [Guid("2F2B0413-1F83-4777-B3B4-38DE3C32DC6B"),
   ClassInterface(ClassInterfaceType.None),
   ComVisible(true)]
   [Obsolete("This class will be withdrawn in the next major release, please use the SOFA or NOVAS31 classes instead")]
   public class Ephemeris : IEphemeris
   {

      private const double DTVEL = 0.01;

      //Ephemeris variables
      private string m_Name; //Name of body
      private Body m_Number; //Number of body
      private bool m_bNumberValid;
      private BodyType m_Type; //Type of body
      private bool m_bTypeValid;
      private orbit m_e; //Elements, etc for minor planets/comets, etc.
                         //Private TL As TraceLogger
                         //gplan variables
      private double[,] ss; // VBConversions Note: Initial value cannot be assigned here since it is non-static.  Assignment has been moved to the class constructors.
      private double[,] cc; // VBConversions Note: Initial value cannot be assigned here since it is non-static.  Assignment has been moved to the class constructors.
      private double[] Args; // VBConversions Note: Initial value cannot be assigned here since it is non-static.  Assignment has been moved to the class constructors.
      private double LP_equinox;
      private double NF_arcsec;
      private double Ea_arcsec;
      private double pA_precession;

      /// <summary>
      /// Create a new Ephemeris component and initialise it
      /// </summary>
      /// <remarks></remarks>
      public Ephemeris()
      {
         // VBConversions Note: Non-static class variable initialization is below.  Class variables cannot be initially assigned non-static values in C#.
         ss = new double[NARGS + 1, 32];
         cc = new double[NARGS + 1, 32];
         Args = new double[NARGS + 1];

         //TL = New TraceLogger("", "KeplerEphemeris")
         //TL.Enabled = True
         //TL.LogMessage("New", "Created")
         m_bTypeValid = false;
         m_Name = ""; //Sentinel
         m_Type = null;
         m_e.ptable.lon_tbl = new double[] { 0.0 }; //Initialise orbit arrays
         m_e.ptable.lat_tbl = new double[] { 0.0 };
      }
      /// <summary>
      /// Semi-major axis (AU)
      /// </summary>
      /// <value>Semi-major axis in AU</value>
      /// <returns>Semi-major axis in AU</returns>
      /// <remarks></remarks>
      public double a
      {
         get
         {
            return m_e.a;
         }
         set
         {
            m_e.a = value;
         }
      }

      /// <summary>
      /// The type of solar system body represented by this instance of the ephemeris engine (enum)
      /// </summary>
      /// <value>The type of solar system body represented by this instance of the ephemeris engine (enum)</value>
      /// <returns>0 for major planet, 1 for minot planet and 2 for comet</returns>
      /// <remarks></remarks>
      public BodyType BodyType
      {
         get
         {
            if (!m_bTypeValid) {
               throw (new Exceptions.ValueNotSetException("KEPLER:BodyType BodyType has not been set"));
            }
            return m_Type;
         }
         set
         {
            m_Type = value;
            m_bTypeValid = true;
         }
      }

      /// <summary>
      /// Orbital eccentricity
      /// </summary>
      /// <value>Orbital eccentricity </value>
      /// <returns>Orbital eccentricity </returns>
      /// <remarks></remarks>
      public double e
      {
         get
         {
            return m_e.ecc;
         }
         set
         {
            m_e.ecc = value;
         }
      }

      /// <summary>
      /// Epoch of osculation of the orbital elements (terrestrial Julian date)
      /// </summary>
      /// <value>Epoch of osculation of the orbital elements</value>
      /// <returns>Terrestrial Julian date</returns>
      /// <remarks></remarks>
      public double Epoch
      {
         get
         {
            return m_e.epoch;
         }
         set
         {
            m_e.epoch = value;
         }
      }

      /// <summary>
      /// Slope parameter for magnitude
      /// </summary>
      /// <value>Slope parameter for magnitude</value>
      /// <returns>Slope parameter for magnitude</returns>
      /// <remarks></remarks>
      public double G
      {
         get
         {
            throw (new Exceptions.ValueNotAvailableException("Kepler:G Read - Magnitude slope parameter calculation not implemented"));
         }
         set
         {
            throw (new Exceptions.ValueNotAvailableException("Kepler:G Write - Magnitude slope parameter calculation not implemented"));
         }
      }

      /// <summary>
      /// Compute rectangular (x/y/z) heliocentric J2000 equatorial coordinates of position (AU) and
      /// velocity (KM/sec.).
      /// </summary>
      /// <param name="tjd">Terrestrial Julian date/time for which position and velocity is to be computed</param>
      /// <returns>Array of 6 values containing rectangular (x/y/z) heliocentric J2000 equatorial
      /// coordinates of position (AU) and velocity (KM/sec.) for the body.</returns>
      /// <remarks>The TJD parameter is the date/time as a Terrestrial Time Julian date. See below for
      /// more info. If you are using ACP, there are functions available to convert between UTC and
      /// Terrestrial time, and for estimating the current value of delta-T. See the Overview page for
      /// the Kepler.Ephemeris class for more information on time keeping systems.</remarks>
      public double[] GetPositionAndVelocity(double tjd)
      {
         double[] posvec = new double[6];
         int[] ai = new int[2];
         double[,] pos = new double[4, 4];
         orbit op = new orbit();
         int i = 0;
         double qjd = 0;
         double[] p = new double[3];

         if (!m_bTypeValid) {
            throw (new Exceptions.ValueNotSetException("Kepler:GetPositionAndVelocity Body type has not been set"));
         }
         //TL.LogMessage("GetPosAndVel", m_Number.ToString)
         if (m_Type == BodyType.MajorPlanet) //MAJOR PLANETS [unimpl. SUN, MOON]
         {
            if (m_Number == Body.Mercury) {
               op = mercury;
            }
            else if (m_Number == Body.Venus) {
               op = venus;
            }
            else if (m_Number == Body.Earth) {
               op = earthplanet;
            }
            else if (m_Number == Body.Mars) {
               op = mars;
            }
            else if (m_Number == Body.Jupiter) {
               op = jupiter;
            }
            else if (m_Number == Body.Saturn) {
               op = saturn;
            }
            else if (m_Number == Body.Uranus) {
               op = uranus;
            }
            else if (m_Number == Body.Neptune) {
               op = neptune;
            }
            else if (m_Number == Body.Pluto) {
               op = pluto;
            }
            else {
               throw (new Utilities.Exceptions.InvalidValueException("Kepler:GetPositionAndVelocity Invalid value for planet number: " + System.Convert.ToString(m_Number)));
            }
         } //MINOR PLANET
         else if (m_Type == BodyType.MinorPlanet) {
            ////TODO: Check elements
            op = m_e;
         } //COMET
         else if (m_Type == BodyType.Comet) {
            ////TODO: Check elements
            op = m_e;
         }
         for (i = 0; i <= 2; i++) {
            qjd = tjd + (i - 1) * DTVEL;
            //TL.LogMessage("GetPosVel", "Before KepCalc")
            KeplerCalc(qjd, op, p);
            //TL.LogMessage("GetPosVel", "After KepCalc")
            pos[i, 0] = p[0];
            pos[i, 1] = p[1];
            pos[i, 2] = p[2];
         }

         //pos(1,x) contains the pos vector
         //pos(0,x) and pos(2,x) are used to determine the velocity based on position change with time!
         for (i = 0; i <= 2; i++) {
            posvec[i] = pos[1, i];
            posvec[3 + i] = System.Convert.ToInt32((pos[2, i] - pos[0, i]) / (2.0 * DTVEL));
         }

         return posvec;
      }

      /// <summary>
      /// Absolute visual magnitude
      /// </summary>
      /// <value>Absolute visual magnitude</value>
      /// <returns>Absolute visual magnitude</returns>
      /// <remarks></remarks>
      public double H
      {
         get
         {
            throw (new Exceptions.ValueNotAvailableException("Kepler:H Read - Visual magnitude calculation not implemented"));
         }
         set
         {
            throw (new Exceptions.ValueNotAvailableException("Kepler:H Write - Visual magnitude calculation not implemented"));
         }
      }

      /// <summary>
      /// The J2000.0 inclination (deg.)
      /// </summary>
      /// <value>The J2000.0 inclination</value>
      /// <returns>Degrees</returns>
      /// <remarks></remarks>
      public double Incl
      {
         get
         {
            return m_e.i;
         }
         set
         {
            m_e.i = value;
         }
      }

      /// <summary>
      /// Mean anomaly at the epoch
      /// </summary>
      /// <value>Mean anomaly at the epoch</value>
      /// <returns>Mean anomaly at the epoch</returns>
      /// <remarks></remarks>
      public double M
      {
         get
         {
            return m_e.M;
         }
         set
         {
            m_e.M = value;
         }
      }

      /// <summary>
      /// Mean daily motion (deg/day)
      /// </summary>
      /// <value>Mean daily motion</value>
      /// <returns>Degrees per day</returns>
      /// <remarks></remarks>
      public double n
      {
         get
         {
            return m_e.dm;
         }
         set
         {
            m_e.dm = value;
         }
      }

      /// <summary>
      /// The name of the body.
      /// </summary>
      /// <value>The name of the body or packed MPC designation</value>
      /// <returns>The name of the body or packed MPC designation</returns>
      /// <remarks>If this instance represents an unnumbered minor planet, Ephemeris.Name must be the
      /// packed MPC designation. For other types, this is for display only.</remarks>
      public string Name
      {
         get
         {
            if (string.IsNullOrEmpty(m_Name)) {
               throw (new Exceptions.ValueNotSetException("KEPLER:Name Name has not been set"));
            }
            return m_Name;
         }
         set
         {
            m_Name = value;
         }
      }

      /// <summary>
      /// The J2000.0 longitude of the ascending node (deg.)
      /// </summary>
      /// <value>The J2000.0 longitude of the ascending node</value>
      /// <returns>Degrees</returns>
      /// <remarks></remarks>
      public double Node
      {
         get
         {
            return m_e.W;
         }
         set
         {
            m_e.W = value;
         }
      }

      /// <summary>
      /// The major or minor planet number
      /// </summary>
      /// <value>The major or minor planet number</value>
      /// <returns>Number or zero if not numbered</returns>
      /// <remarks></remarks>
      public Body Number
      {
         get
         {
            if (!m_bNumberValid) {
               throw (new Exceptions.ValueNotSetException("KEPLER:Number Planet number has not been set"));
            }
            return m_Number;
         }
         set
         {
            m_Number = value;
         }
      }

      /// <summary>
      /// Orbital period (years)
      /// </summary>
      /// <value>Orbital period</value>
      /// <returns>Years</returns>
      /// <remarks></remarks>
      public double P
      {
         get
         {
            throw (new Exceptions.ValueNotAvailableException("Kepler:P Read - Orbital period calculation not implemented"));
         }
         set
         {
            throw (new Exceptions.ValueNotAvailableException("Kepler:P Write - Orbital period calculation not implemented"));
         }
      }

      /// <summary>
      /// The J2000.0 argument of perihelion (deg.)
      /// </summary>
      /// <value>The J2000.0 argument of perihelion</value>
      /// <returns>Degrees</returns>
      /// <remarks></remarks>
      public double Peri
      {
         get
         {
            return m_e.wp;
         }
         set
         {
            m_e.wp = value;
         }
      }

      /// <summary>
      /// Perihelion distance (AU)
      /// </summary>
      /// <value>Perihelion distance</value>
      /// <returns>AU</returns>
      /// <remarks></remarks>
      public double q
      {
         get
         {
            return m_e.a;
         }
         set
         {
            m_e.a = value;
         }
      }

      /// <summary>
      /// Reciprocal semi-major axis (1/AU)
      /// </summary>
      /// <value>Reciprocal semi-major axis</value>
      /// <returns>1/AU</returns>
      /// <remarks></remarks>
      public double z
      {
         get
         {
            return 1.0 / m_e.a;
         }
         set
         {
            m_e.a = 1.0 / value;
         }
      }
   }
}