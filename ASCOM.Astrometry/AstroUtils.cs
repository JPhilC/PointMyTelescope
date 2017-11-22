

namespace ASCOM.Astrometry
{
   /// <summary>
   /// Class providing a suite of tested astronomy support functions to save develpment effort and provide consistant behaviour.
   /// </summary>
   /// <remarks>
   /// A number of these routines are provided to support migration from the Astro32.dll. Unlike Astro32, these routines will work in
   /// both 32bit and 64bit applications.
   /// </remarks>
   public class AstroUtils 
   {


      public double DeltaUT(double JulianDate)
      {
         double DUT1 = 0;

         DUT1 = System.Convert.ToDouble(GlobalItems.TAI_UTC_OFFSET + GlobalItems.TT_TAI_OFFSET - DeltatCode.DeltaTCalc(JulianDate));
         return DUT1;
      }



   }
}