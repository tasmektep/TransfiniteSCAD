using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SCAD
{


    public enum Blending_Methods
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

    /// <summary>
    ///  Input  : 
    ///  Output :
    /// </summary>
    class BlendingFunctions
    {






       


        public BlendingFunctions(int distance, Blending_Methods method)
        {
            if( method == Blending_Methods.Side_Blending)
            {

            }

            if (method == Blending_Methods.Special_Side_Blending)
            {

            }

            if (method == Blending_Methods.Corner_Blending)
            {

            }
        }

        //Special Side Blending


        //Corner Blending


        //Side Blending (Kato)

    }
}
