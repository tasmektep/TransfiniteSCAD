using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SCAD
{
    /// <summary>
    ///  Input  : 
    ///  Output :
    /// </summary>
    class BlendingFunctions
    {






        public enum Method
        {
            //
            // Summary:
            //      Uses Side Blending functions
            Side_Blending = 0,
            //
            // Summary:
            //      Uses Corner Blending functions
            Corner_Blending,
            //
            // Summary:
            //      Uses Special Side Blending functions
            Special_Side_Blending


        }


        public BlendingFunctions(int distance, Method method)
        {
            if( method == Method.Side_Blending)
            {

            }

            if (method == Method.Special_Side_Blending)
            {

            }

            if (method == Method.Corner_Blending)
            {

            }
        }

        //Special Side Blending


        //Corner Blending


        //Side Blending (Kato)

    }
}
