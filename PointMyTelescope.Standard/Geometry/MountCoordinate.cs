﻿using ASCOM.Astrometry.Transform;
// using ASCOM.Utilities;

namespace Lunatic.Core.Geometry
{

   /// <summary>
   /// A class to tie together an equatorial coordinate, 
   /// calculated/theoretical mount axis positions at a given time
   /// and optionally the actual observered axis positions.
   /// </summary>
   public class MountCoordinate
   {
      public enum MasterCoordinateEnum
      {
         Equatorial,
         AltAzimuth
      }

      private EquatorialCoordinate _Equatorial = new EquatorialCoordinate();
      private AltAzCoordinate _AltAzimuth = new AltAzCoordinate();
      private IntegerPair _AxesPosition = new IntegerPair();
      private double _AxisJulianTimeUTC;
      private double _LocalJulianTimeUTC;

      /// <summary>
      /// Held for reference so that when a refresh is requested we know which coordinate
      /// is the master.
      /// </summary>
      private MasterCoordinateEnum _MasterCoordinate;

      public MasterCoordinateEnum MasterCoordinate
      {
         get
         {
            return _MasterCoordinate;
         }
         private set
         {
            _MasterCoordinate = value;
         }
      }

      public EquatorialCoordinate Equatorial
      {
         get
         {
            return _Equatorial;
         }
         //private set
         //{
         //   _Equatorial = value;
         //}
      }

      public AltAzCoordinate AltAzimuth
      {
         get
         {
            return _AltAzimuth;
         }
         set
         {
            _AltAzimuth = value;
         }
      }


      public IntegerPair ObservedAxes
      {
         get
         {
            return _AxesPosition;
         }
         //private set
         //{
         //   _AxesPosition = value;
         //}
      }

      /// <summary>
      /// The time at which the axis values were determined.
      /// </summary>
      public double AxisJulianTimeUTC
      {
         get
         {
            return _AxisJulianTimeUTC;
         }
         //private set
         //{
         //   _AxisJulianTimeUTC = value;
         //}
      }

      /// <summary>
      /// Julian UTC time of last update/synch
      /// </summary>
      public double LocalJulianTimeUTC
      {
         get
         {
            return _LocalJulianTimeUTC;
         }
         private set
         {
            _LocalJulianTimeUTC = value;
         }
      }


      #region Constructors ...
      /// <summary>
      /// Initialise a mount coordinate with Ra/Dec strings 
      /// </summary>
      /// <param name="ra">A right ascension string</param>
      /// <param name="dec">declination string</param>
      /// <param name="localTime">The local time of the observation</param>
      public MountCoordinate(string ra, string dec) : this(new EquatorialCoordinate(ra, dec))
      {
         _MasterCoordinate = MasterCoordinateEnum.Equatorial;
      }

      /// <summary>
      /// Simple initialisation with an equatorial coordinate
      /// </summary>
      public MountCoordinate(EquatorialCoordinate equatorial)
      {
         _Equatorial = equatorial;
         _MasterCoordinate = MasterCoordinateEnum.Equatorial;
      }

      /// <summary>
      /// Simple initialisation with an altAzimuth coordinate
      /// </summary>
      public MountCoordinate(AltAzCoordinate altAz)
      {
         _AltAzimuth = altAz;
         _MasterCoordinate = MasterCoordinateEnum.AltAzimuth;
      }

      /// <summary>
      /// Initialisation with an equatorial coordinate, a transform instance and the local julian time (corrected)
      /// which then means that the AltAzimunth at the time is available.
      /// </summary>
      public MountCoordinate(EquatorialCoordinate equatorial, Transform transform, double localJulianTimeUTC) : this(equatorial)
      {
         _AltAzimuth = this.GetAltAzimuth(transform, localJulianTimeUTC);
         _LocalJulianTimeUTC = localJulianTimeUTC;
      }

      /// <summary>
      /// Initialise a mount coordinate with Ra/Dec strings and axis positions in radians.
      /// </summary>
      /// <param name="altAz">The AltAzimuth coordinate for the mount</param>
      /// <param name="suggested">The suggested position for the axes (e.g. via a star catalogue lookup)</param>
      /// <param name="localTime">The local time of the observation</param>
      public MountCoordinate(string ra, string dec, IntegerPair axisPosition, Transform transform, double localJulianTimeUTC) : this(new EquatorialCoordinate(ra, dec))
      {
         _Equatorial = new EquatorialCoordinate(ra, dec);
         _AltAzimuth = this.GetAltAzimuth(transform, localJulianTimeUTC);
         _AxesPosition = axisPosition;
         _AxisJulianTimeUTC = localJulianTimeUTC;
         _MasterCoordinate = MasterCoordinateEnum.Equatorial;
      }


      /// <summary>
      /// Initialisation with an equatorial coordinate, a transform instance and the local julian time (corrected)
      /// which then means that the AltAzimunth at the time is available.
      /// </summary>
      public MountCoordinate(AltAzCoordinate altAz, Transform transform, double localJulianTimeUTC) : this(altAz)
      {
         _LocalJulianTimeUTC = localJulianTimeUTC;
         _Equatorial = this.GetEquatorial(transform, localJulianTimeUTC);
      }

      #endregion

      /// <summary>
      /// Returns the AltAzimuth coordinate for the equatorial using the values
      /// currently set in the passed Transform instance.
      /// </summary>
      /// <param name="transform"></param>
      /// <returns></returns>
      public AltAzCoordinate GetAltAzimuth(Transform transform)
      {
         transform.SetTopocentric(_Equatorial.RightAscension, _Equatorial.Declination);
         transform.Refresh();
         AltAzCoordinate coord = new AltAzCoordinate(transform.ElevationTopocentric, transform.AzimuthTopocentric);
         return coord;
      }

      /// <summary>
      /// Returns the AltAzimuth coordinate for the equatorial using the values
      /// currently set in the passed Transform instance.
      /// </summary>
      /// <param name="transform"></param>
      /// <returns></returns>
      public AltAzCoordinate GetAltAzimuth(Transform transform, double julianDateUTC)
      {
         transform.JulianDateUTC = julianDateUTC;
         transform.SetTopocentric(_Equatorial.RightAscension, _Equatorial.Declination);
         transform.Refresh();
         AltAzCoordinate coord = new AltAzCoordinate(transform.ElevationTopocentric, transform.AzimuthTopocentric);
         return coord;
      }

      private void RefreshAltAzimuth(Transform transform, double julianDateUTC)
      {
         transform.JulianDateUTC = julianDateUTC;
         transform.SetTopocentric(_Equatorial.RightAscension, _Equatorial.Declination);
         transform.Refresh();
         _AltAzimuth.Altitude.Value = transform.ElevationTopocentric;
         _AltAzimuth.Azimuth.Value = transform.AzimuthTopocentric;
      }

      /// <summary>
      /// Returns the RADec coordinate for the observed AltAzimuth using the values
      /// currently set in the passed Transform instance.
      /// </summary>
      /// <param name="transform"></param>
      /// <returns></returns>
      public EquatorialCoordinate GetEquatorial(Transform transform, double julianDateUTC)
      {
         transform.JulianDateUTC = julianDateUTC;
         transform.SetAzimuthElevation(_AltAzimuth.Azimuth, _AltAzimuth.Altitude);
         transform.Refresh();
         EquatorialCoordinate coord = new EquatorialCoordinate(transform.RATopocentric, transform.DECTopocentric);
         return coord;
      }


      public void SetObservedAxis(IntegerPair axisPosition, double observationTime)
      {
         _AxesPosition = axisPosition;
         _AxisJulianTimeUTC = observationTime;
      }

      public void Refresh(Transform transform, double julianDateUTC)
      {
         transform.JulianDateUTC = julianDateUTC;
         if (_MasterCoordinate == MasterCoordinateEnum.Equatorial) {
            // Update the AltAzimuth
            transform.SetTopocentric(_Equatorial.RightAscension, _Equatorial.Declination);
            transform.Refresh();
            _AltAzimuth.Altitude.Value = transform.ElevationTopocentric;
            _AltAzimuth.Azimuth.Value = transform.AzimuthTopocentric;
         }
         else {
            // Update the Equatorial
            transform.SetAzimuthElevation(_AltAzimuth.Azimuth, _AltAzimuth.Altitude);
            transform.Refresh();
            _Equatorial.RightAscension.Value = transform.RATopocentric;
            _Equatorial.Declination.Value = transform.DECTopocentric;
         }
      }

      public void Refresh(EquatorialCoordinate equatorial, IntegerPair axisPosition, Transform transform, double localJulianTimeUTC)
      {
         _Equatorial.RightAscension.Value = equatorial.RightAscension.Value;
         _Equatorial.Declination.Value = equatorial.Declination.Value;
         _AxesPosition = axisPosition;

         _LocalJulianTimeUTC = localJulianTimeUTC;
         RefreshAltAzimuth(transform, localJulianTimeUTC);
         _MasterCoordinate = MasterCoordinateEnum.Equatorial;

      }


   }

}
