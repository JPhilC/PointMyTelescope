using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace Lunatic.Core.Geometry
{

   /// <summary>
   /// A structure to represent telecope mount axis positions
   /// </summary>
   public class AxisPosition : DataObjectBase
   {
      private const string RAAXIS_HAS_ERRORS = "RA axis has errors.";
      private const string RAAXIS_RANGE_ERROR = "RA axis value be greater than or equal to 0 and less than 360 degrees.";
      private const string DECAXIS_HAS_ERRORS = "Dec axis has errors.";
      private const string DECAXIS_RANGE_ERROR = "Dec axis value be greater than or equal to 0 and less than 360 degrees.";

      private Angle _RAAxis;
      private Angle _DecAxis;

      public Angle RAAxis
      {
         get
         {
            return _RAAxis;
         }
         private set
         {
            if (_RAAxis != null) {
               WeakEventManager<Angle, DataErrorsChangedEventArgs>.RemoveHandler(_RAAxis, "ErrorsChanged", RAAxis_ErrorsChanged);
               WeakEventManager<Angle, PropertyChangedEventArgs>.RemoveHandler(_RAAxis, "PropertyChanged", RAAxis_PropertyChanged);
            }
            _RAAxis = value;
            _RAAxis.Format = AngularFormat.DecimalDegrees;
            if (_RAAxis != null) {
               WeakEventManager<Angle, DataErrorsChangedEventArgs>.AddHandler(_RAAxis, "ErrorsChanged", RAAxis_ErrorsChanged);
               WeakEventManager<Angle, PropertyChangedEventArgs>.AddHandler(_RAAxis, "PropertyChanged", RAAxis_PropertyChanged);
            }
            RaisePropertyChanged("RAAxis");
         }
      }

      private void RAAxis_PropertyChanged(object sender, PropertyChangedEventArgs e)
      {
         if (e.PropertyName == "Value") {
            if (RAAxis.Value < 0.0 || RAAxis >= 360.0) {
               AddError("RAAxis", RAAXIS_RANGE_ERROR);
            }
            else {
               RemoveError("RAAxis", RAAXIS_RANGE_ERROR);
            }
         }
         RaisePropertyChanged("RAAxis");
      }

      private void RAAxis_ErrorsChanged(object sender, DataErrorsChangedEventArgs e)
      {
         if (RAAxis.HasErrors) {
            AddError("RAAxis", RAAXIS_HAS_ERRORS);
         }
         else {
            RemoveError("RAAxis", RAAXIS_HAS_ERRORS);
         }
      }


      public Angle DecAxis
      {
         get
         {
            return _DecAxis;
         }
         private set
         {
            if (_DecAxis != null) {
               WeakEventManager<Angle, DataErrorsChangedEventArgs>.RemoveHandler(_DecAxis, "ErrorsChanged", DecAxis_ErrorsChanged);
               WeakEventManager<Angle, PropertyChangedEventArgs>.RemoveHandler(_DecAxis, "PropertyChanged", DecAxis_PropertyChanged);
            }
            _DecAxis = value;
            _DecAxis.Format = AngularFormat.DecimalDegrees;
            if (_DecAxis != null) {
               WeakEventManager<Angle, DataErrorsChangedEventArgs>.AddHandler(_DecAxis, "ErrorsChanged", DecAxis_ErrorsChanged);
               WeakEventManager<Angle, PropertyChangedEventArgs>.AddHandler(_DecAxis, "PropertyChanged", DecAxis_PropertyChanged);
            }
            RaisePropertyChanged("DecAxis");
         }
      }
      private void DecAxis_PropertyChanged(object sender, PropertyChangedEventArgs e)
      {
         if (e.PropertyName == "Value") {
            if (DecAxis.Value < 0.0 || DecAxis >= 360.0) {
               AddError("DecAxis", DECAXIS_RANGE_ERROR);
            }
            else {
               RemoveError("DecAxis", DECAXIS_RANGE_ERROR);
            }
         }
         RaisePropertyChanged("DecAxis");
      }

      private void DecAxis_ErrorsChanged(object sender, DataErrorsChangedEventArgs e)
      {
         if (DecAxis.HasErrors) {
            AddError("DecAxis", DECAXIS_HAS_ERRORS);
         }
         else {
            RemoveError("DecAxis", DECAXIS_HAS_ERRORS);
         }
      }

      public int AxisCount
      {
         get
         {
            return 2;
         }
      }

      public AxisPosition()
      {
         RAAxis = new Angle(0.0);
         DecAxis = new Angle(0.0);
      }

      /// <summary>
      /// Initialise the Axis positions
      /// </summary>
      /// <param name="raPosition">RA Axis position in degrees</param>
      /// <param name="decPosition">Dec Axis position in degrees</param>
      public AxisPosition(string raPosition, string decPosition)
      {
         RAAxis = new Angle(raPosition);
         DecAxis = new Angle(decPosition);
      }
      public AxisPosition(double raRadians, double decRadians) : this()
      {
         RAAxis.Radians = raRadians;
         DecAxis.Radians = decRadians;
      }

      public void Update(string axisPositions)
      {
         System.Diagnostics.Debug.WriteLine("Axis positions: " + axisPositions);
         string[] positions = axisPositions.Split(',');
         try {
            RAAxis.Value = double.Parse(positions[0]);
            DecAxis.Value = double.Parse(positions[1]);
         }
         catch {
            throw new ArgumentException("Badly formed axis position string");
         }
      }


      /// <summary>
      /// Returns the index axis value in Radians.
      /// </summary>
      /// <param name="index"></param>
      /// <returns></returns>
      public double this[int index]
      {
         get
         {
            if (index < 0 || index > 1) {
               throw new ArgumentOutOfRangeException();
            }
            return (index == 0 ? RAAxis.Radians : DecAxis.Radians);
         }
         set
         {
            if (index < 0 || index > 1) {
               throw new ArgumentOutOfRangeException();
            }
            if (index == 0) {
               RAAxis.Radians = value;
            }
            else {
               DecAxis.Radians = value;
            }
         }
      }

      /// <summary>
      /// Compares the two specified sets of Axis positions.
      /// </summary>
      public static bool operator ==(AxisPosition pos1, AxisPosition pos2)
      {
         return (pos1.RAAxis.Radians == pos2.RAAxis.Radians && pos1.DecAxis.Radians == pos2.DecAxis.Radians);
      }

      public static bool operator !=(AxisPosition pos1, AxisPosition pos2)
      {
         return !(pos1 == pos2);
      }

      public static AxisPosition operator -(AxisPosition pos1, AxisPosition pos2)
      {
         return new AxisPosition(pos1.RAAxis.Radians - pos2.RAAxis.Radians, pos1.DecAxis.Radians - pos2.DecAxis.Radians);
      }

      public static AxisPosition operator +(AxisPosition pos1, AxisPosition pos2)
      {
         return new AxisPosition(pos1.RAAxis.Radians + pos2.RAAxis.Radians, pos1.DecAxis.Radians + pos2.DecAxis.Radians);
      }

      public override int GetHashCode()
      {
         unchecked // Overflow is fine, just wrap
         {
            int hash = 17;
            // Suitable nullity checks etc, of course :)
            hash = hash * 23 + RAAxis.GetHashCode();
            hash = hash * 23 + DecAxis.GetHashCode();
            return hash;
         }
      }

      public override bool Equals(object obj)
      {
         return (obj is AxisPosition
                 && this == (AxisPosition)obj);
      }


      public bool Equals(AxisPosition obj, double toleranceRadians)
      {
         return ((Math.Abs(obj.RAAxis.Radians - this.RAAxis.Radians) < toleranceRadians)
            && (Math.Abs(obj.DecAxis.Radians - this.DecAxis.Radians) < toleranceRadians));
      }

      public override string ToString()
      {
         return string.Format("RAAxis = {0} Radians, DecAxis = {1} Radians", RAAxis.Radians, DecAxis.Radians);
      }
      public string ToDegreesString()
      {
         return string.Format("{0},{1}", RAAxis.Value, DecAxis.Value);
      }
      public string ToRadiansString()
      {
         return string.Format("{0},{1}", RAAxis.Radians, DecAxis.Radians);
      }

   }

}
