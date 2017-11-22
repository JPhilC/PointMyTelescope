using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.ComponentModel;
using Newtonsoft.Json;

namespace Lunatic.Core.Geometry
{
   /// <summary>
   /// A structure to represent an EquatorialCoordinate
   /// </summary>
   public class EquatorialCoordinate : DataObjectBase, IEquatable<EquatorialCoordinate>
   {
      private HourAngle _RightAscension;
      private Angle _Declination;
      //private DateTime _ObservedWhen;
      //private Angle _Longitude;
      
      private const string RIGHTASCENSION_HAS_ERRORS = "Right Ascension errors exist.";
      private const string RIGHTASCENSION_RANGE_ERROR = "Right Ascension must be greater than or equal to 0.0 and less than 24.0.";
      private const string DECLINATION_HAS_ERRORS = "Declination errors exist.";
      private const string DECLINATION_RANGE_ERROR = "Declination must lie between -90.0 and 90.0";

      [JsonProperty]
      public HourAngle RightAscension
      {
         get
         {
            return _RightAscension;
         }
         private set
         {
            if (_RightAscension != null) {
               WeakEventManager<HourAngle, DataErrorsChangedEventArgs>.RemoveHandler(_RightAscension, "ErrorsChanged", RightAscension_ErrorsChanged);
               WeakEventManager<HourAngle, PropertyChangedEventArgs>.RemoveHandler(_RightAscension, "PropertyChanged", RightAscension_PropertyChanged);
            }
            _RightAscension = value;
            if (_RightAscension != null) {
               WeakEventManager<HourAngle, DataErrorsChangedEventArgs>.AddHandler(_RightAscension, "ErrorsChanged", RightAscension_ErrorsChanged);
               WeakEventManager<HourAngle, PropertyChangedEventArgs>.AddHandler(_RightAscension, "PropertyChanged", RightAscension_PropertyChanged);
            }
            RaisePropertyChanged("RightAscension");
         }
      }

      private void RightAscension_PropertyChanged(object sender, PropertyChangedEventArgs e)
      {
         if (e.PropertyName == "Value") {
            if (RightAscension.Value < 0.0 || RightAscension >= 24.0) {
               AddError("RightAscension", RIGHTASCENSION_RANGE_ERROR);
            }
            else {
               RemoveError("RightAscension", RIGHTASCENSION_RANGE_ERROR);
            }
         }
         RaisePropertyChanged("RightAscension");
      }

      private void RightAscension_ErrorsChanged(object sender, DataErrorsChangedEventArgs e)
      {
         if (RightAscension.HasErrors) {
            AddError("RightAscension", RIGHTASCENSION_HAS_ERRORS);
         }
         else {
            RemoveError("RightAscension", RIGHTASCENSION_HAS_ERRORS);
         }
      }

      public Angle Declination
      {
         get
         {
            return _Declination;
         }
         private set
         {
            if (_Declination != null) {
               WeakEventManager<Angle, DataErrorsChangedEventArgs>.RemoveHandler(_Declination, "ErrorsChanged", Declination_ErrorsChanged);
               WeakEventManager<Angle, PropertyChangedEventArgs>.RemoveHandler(_Declination, "PropertyChanged", Declination_PropertyChanged);
            }
            _Declination = value;
            if (_Declination != null) {
               WeakEventManager<Angle, DataErrorsChangedEventArgs>.AddHandler(_Declination, "ErrorsChanged", Declination_ErrorsChanged);
               WeakEventManager<Angle, PropertyChangedEventArgs>.AddHandler(_Declination, "PropertyChanged", Declination_PropertyChanged);
            }
            RaisePropertyChanged("Declination");
         }
      }

      private void Declination_PropertyChanged(object sender, PropertyChangedEventArgs e)
      {
         if (e.PropertyName == "Value") {
            if (Declination.Value < -90.0 || Declination >= 90.0) {
               AddError("Declination", DECLINATION_RANGE_ERROR);
            }
            else {
               RemoveError("Declination", DECLINATION_RANGE_ERROR);
            }
         }
         RaisePropertyChanged("Declination");
      }

      private void Declination_ErrorsChanged(object sender, DataErrorsChangedEventArgs e)
      {
         if (Declination.HasErrors) {
            AddError("Declination", DECLINATION_HAS_ERRORS);
         }
         else {
            RemoveError("Declination", DECLINATION_HAS_ERRORS);
         }
      }

      public EquatorialCoordinate()   // , double longitude, DateTime observedTime)
      {
         RightAscension = new HourAngle(0.0);
         Declination = new Angle(0.0);
      }

      public EquatorialCoordinate(double rightAscension, double declination) : this()   // , double longitude, DateTime observedTime)
      {
         if (rightAscension < 0 || rightAscension > 24.0) { throw new ArgumentOutOfRangeException(RIGHTASCENSION_RANGE_ERROR); }
         if (declination < -90 || declination > 90) { throw new ArgumentOutOfRangeException(DECLINATION_RANGE_ERROR); }
         RightAscension.Value = rightAscension;
         Declination.Value = declination;
      }


      public EquatorialCoordinate(HourAngle rightAscension, Angle declination)    // , Angle longitude, DateTime observedTime)
      {
         if (rightAscension.Value < 0 || rightAscension.Value > 24.0) { throw new ArgumentOutOfRangeException(RIGHTASCENSION_RANGE_ERROR); }
         if (declination.Value < -90 || declination.Value > 90) { throw new ArgumentOutOfRangeException(DECLINATION_RANGE_ERROR); }
         RightAscension = rightAscension;
         Declination = declination;
      }

      #region Operator overloads ...
      /// <summary>
      /// Compares the two specified sets of Axis positions.
      /// </summary>
      public static bool operator ==(EquatorialCoordinate pos1, EquatorialCoordinate pos2)
      {
         if (ReferenceEquals(pos1, pos2)) {
            return true;
         }

         if (ReferenceEquals(pos1, null)) {
            return false;
         }
         if (ReferenceEquals(pos2, null)) {
            return false;
         }
         return (pos1.RightAscension.Value == pos2.RightAscension.Value && pos1.Declination.Value == pos2.Declination.Value);
      }

      public static bool operator !=(EquatorialCoordinate pos1, EquatorialCoordinate pos2)
      {
         return !(pos1 == pos2);
      }

      public override int GetHashCode()
      {
         unchecked // Overflow is fine, just wrap
         {
            int hash = 17;
            // Suitable nullity checks etc, of course :)
            hash = hash * 23 + RightAscension.GetHashCode();
            hash = hash * 23 + Declination.GetHashCode();
            return hash;
         }
      }

      public bool Equals(EquatorialCoordinate other)
      {
         if (other == null)
            return false;

         return (this.RightAscension.Value == other.RightAscension.Value && this.Declination.Value == other.Declination.Value);
      }

      public override bool Equals(object obj)
      {
         if (ReferenceEquals(null, obj)) {
            return false;
         }
         if (ReferenceEquals(this, obj)) {
            return true;
         }

         return obj.GetType() == GetType() && Equals((EquatorialCoordinate)obj);
      }


      public static EquatorialCoordinate operator -(EquatorialCoordinate pos1, EquatorialCoordinate pos2)
      {
         return new EquatorialCoordinate(pos1.RightAscension - pos2.RightAscension, pos1.Declination - pos2.Declination);
      }

      public static EquatorialCoordinate operator +(EquatorialCoordinate pos1, EquatorialCoordinate pos2)
      {
         return new EquatorialCoordinate(pos1.RightAscension + pos2.RightAscension, pos1.Declination + pos2.Declination);
      }


      public override string ToString()
      {
         return string.Format("{0}/{1}", _RightAscension, _Declination);
      }
      #endregion

      public CarteseanCoordinate ToCartesean(Angle latitude, bool affineTaki = true)
      {
         CarteseanCoordinate cartCoord;
         if (affineTaki) {
            // Get Polar (or should than be get AltAzimuth) from Equatorial coordinate (formerly call to EQ_SphericalPolar)
            AltAzCoordinate polar = AstroConvert.GetAltAz(this, latitude);
            // Get  Cartesean from Polar (formerly call to EQ_Polar2Cartes)
            cartCoord = polar.ToCartesean();
         }
         else {
            cartCoord = new CarteseanCoordinate(this.RightAscension.Radians, this.Declination.Radians, 1.0);
         }
         return cartCoord;
      }

      
   }

}
