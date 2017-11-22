using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace Lunatic.Core
{

   /// <summary>
   /// A structure to represent a pair of doubles related to RA and Dec
   /// </summary>
   public struct DoublePair
   {

      public double RA { get; set; }
      public double Dec { get; set; }

      /// <summary>
      /// Initialise the Encoder positions
      /// </summary>
      /// <param name="raPosition">RA Encoder position in degrees</param>
      /// <param name="decPosition">Dec Encoder position in degrees</param>
      public DoublePair(string raPosition, string decPosition)
      {
         RA = double.Parse(raPosition);
         Dec = double.Parse(decPosition);
      }
      public DoublePair(double raEncoder, double decEncoder)
      {
         RA = raEncoder;
         Dec = decEncoder;
      }

      /// <summary>
      /// Returns the index encoder value.
      /// </summary>
      /// <param name="index"></param>
      /// <returns></returns>
      public double this[double index]
      {
         get
         {
            if (index < 0 || index > 1) {
               throw new ArgumentOutOfRangeException();
            }
            return (index == 0 ? RA : Dec);
         }
         set
         {
            if (index < 0 || index > 1) {
               throw new ArgumentOutOfRangeException();
            }
            if (index == 0) {
               RA = value;
            }
            else {
               Dec = value;
            }
         }
      }

      /// <summary>
      /// Compares the two specified sets of Encoder positions.
      /// </summary>
      public static bool operator ==(DoublePair pos1, DoublePair pos2)
      {
         return (pos1.RA == pos2.RA && pos1.Dec == pos2.Dec);
      }

      public static bool operator !=(DoublePair pos1, DoublePair pos2)
      {
         return !(pos1 == pos2);
      }

      public static DoublePair operator -(DoublePair pos1, DoublePair pos2)
      {
         return new DoublePair(pos1.RA - pos2.RA, pos1.Dec - pos2.Dec);
      }

      public static DoublePair operator +(DoublePair pos1, DoublePair pos2)
      {
         return new DoublePair(pos1.RA + pos2.RA, pos1.Dec + pos2.Dec);
      }

      public static DoublePair operator *(DoublePair pos1, DoublePair pos2)
      {
         return new DoublePair(pos1.RA * pos2.RA, pos1.Dec * pos2.Dec);
      }

      public static DoublePair operator /(DoublePair pos1, DoublePair pos2)
      {
         return new DoublePair(pos1.RA / pos2.RA, pos1.Dec / pos2.Dec);
      }

      public static DoublePair operator /(DoublePair pos1, double divisor)
      {
         return new DoublePair(pos1.RA / divisor, pos1.Dec / divisor);
      }


      // Implicit cast override
      public static implicit operator DoublePair(IntegerPair i)
      {
         return new DoublePair(i.RA, i.Dec);
      }

      public IntegerPair ToIntegerPair()
      {
         return new IntegerPair((int)Math.Round(this.RA, 0, MidpointRounding.AwayFromZero),
            (int)Math.Round(this.Dec, 0, MidpointRounding.AwayFromZero));
      }

      public DoublePair Range2Pi()
      {
         return new DoublePair(AstroConvert.Range2Pi(this.RA), AstroConvert.Range2Pi(this.Dec));
      }

      public override int GetHashCode()
      {
         unchecked // Overflow is fine, just wrap
         {
            int hash = 17;
            // Suitable nullity checks etc, of course :)
            hash = hash * 23 + RA.GetHashCode();
            hash = hash * 23 + Dec.GetHashCode();
            return hash;
         }
      }

      public override bool Equals(object obj)
      {
         return (obj is DoublePair
                 && this == (DoublePair)obj);
      }


      public bool Equals(DoublePair obj, double tolerance)
      {
         return ((Math.Abs(obj.RA - this.RA) < tolerance)
            && (Math.Abs(obj.Dec - this.Dec) < tolerance));
      }

      public override string ToString()
      {
         return string.Format("RA = {0}, Dec = {1}", RA, Dec);
      }

      public string ToInterfaceString()
      {
         return string.Format("{0},{1}", RA, Dec);
      }
   }

}
