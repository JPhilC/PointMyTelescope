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
   /// A structure to represent a pair of integers related to RA and Dec.
   /// </summary>
   public struct IntegerPair
   {

      public int RA { get; set; }
      public int Dec { get; set; }

      /// <summary>
      /// Initialise the Encoder positions
      /// </summary>
      /// <param name="raPosition">RA Encoder position in degrees</param>
      /// <param name="decPosition">Dec Encoder position in degrees</param>
      public IntegerPair(string raPosition, string decPosition)
      {
         RA = int.Parse(raPosition);
         Dec = int.Parse(decPosition);
      }

      public IntegerPair(int raEncoder, int decEncoder)
      {
         RA = raEncoder;
         Dec = decEncoder;
      }
      public IntegerPair(double raValue, double decValue)
      {
         RA = (int)Math.Round(raValue, 0, MidpointRounding.AwayFromZero);
         Dec = (int)Math.Round(decValue, 0, MidpointRounding.AwayFromZero);
      }

      /// <summary>
      /// Returns the index encoder value.
      /// </summary>
      /// <param name="index"></param>
      /// <returns></returns>
      public int this[int index]
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
      public static bool operator ==(IntegerPair pos1, IntegerPair pos2)
      {
         return (pos1.RA == pos2.RA && pos1.Dec == pos2.Dec);
      }

      public static bool operator !=(IntegerPair pos1, IntegerPair pos2)
      {
         return !(pos1 == pos2);
      }

      public static IntegerPair operator -(IntegerPair pos1, IntegerPair pos2)
      {
         return new IntegerPair(pos1.RA - pos2.RA, pos1.Dec - pos2.Dec);
      }

      public static IntegerPair operator +(IntegerPair pos1, IntegerPair pos2)
      {
         return new IntegerPair(pos1.RA + pos2.RA, pos1.Dec + pos2.Dec);
      }

      public static IntegerPair operator *(IntegerPair pos1, IntegerPair pos2)
      {
         return new IntegerPair(pos1.RA * pos2.RA, pos1.Dec * pos2.Dec);
      }

      public static IntegerPair operator /(IntegerPair pos1, double divisor)
      {
         return new IntegerPair((int)Math.Round(pos1.RA/divisor, 0, MidpointRounding.AwayFromZero), 
                                 (int)Math.Round(pos1.Dec/divisor, 0, MidpointRounding.AwayFromZero));
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
         return (obj is IntegerPair
                 && this == (IntegerPair)obj);
      }


      public bool Equals(IntegerPair obj, double toleranceRadians)
      {
         return ((Math.Abs(obj.RA - this.RA) < toleranceRadians)
            && (Math.Abs(obj.Dec - this.Dec) < toleranceRadians));
      }

      public IntegerPair AbsOffsetTo(IntegerPair other)
      {
         return new IntegerPair(Math.Abs(this.RA - other.RA), Math.Abs(this.Dec - other.Dec));
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
