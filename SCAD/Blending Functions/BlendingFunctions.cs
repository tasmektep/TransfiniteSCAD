using Grasshopper.Kernel.Geometry;
using Rhino.Input.Custom;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Media3D;
namespace SCAD
{
    public enum Blending_Method
    {
        //
        // Summary:
        //      Uses Side Blending functions
        Side_Blending = 0,
        //
        // Summary:
        //      Uses Corner Blending functions
        Corner_Blending = 1,
        //
        // Summary:
        //      Uses Special Side Blending functions
        Special_Side_Blending = 2,
    }
    /// <summary>
    ///  Input  :
    ///  Output :
    /// </summary>
    class BlendingFunctions
    {
        /// <summary>
        ///
        /// </summary>
        /// <param name="d"></param>
        /// <param name="method"></param>
        /// <returns></returns>
        Blending_Method method = 0;
        //identify employed blending function
        public BlendingFunctions(Blending_Method method)
        {
            this.method = method;
        }
        /// <summary>
        /// For a given distances from edges for a point, blending weight of point for all edges are listed by order
        /// </summary>
        /// <param name="d"></param> distances
        /// <returns></returns> ordered blending weights
        public List<double> GetBlending(List<double> d)
        {
            List<double> Mus = new List<double>();
            if (method == Blending_Method.Side_Blending)
            {
                for (int i = 0; i < d.Count; i++)
                {
                    double Kapa = 0;
                    double Numerator = 0;
                    double denominator = 0;
                    List<int> Numerator1 = new List<int>();
                    List<int> Numerator2 = new List<int>();
                    Numerator1.Add(i);
                    Numerator2.Add(i);
                    //If boundary curve ID equals to 0,its previous neighbour will be the last curve
                    if (i == 0)
                        Numerator1.Add(d.Count - 1);
                    else
                        Numerator1.Add(i - 1);
                    //If boundary curve ID is the last curve,its next neighbour will be first (0) curve
                    if (i == d.Count - 1)
                        Numerator1.Add(0);
                    else
                        Numerator1.Add(i + 1);
                    Numerator = ProductFunction(d, Numerator1) + ProductFunction(d, Numerator2);
                    for (int j = 0; j < d.Count; j++)
                    {
                        List<int> Denominator1 = new List<int>();
                        Denominator1.Add(j);
                        //If boundary curve ID equals to 0,its previous neighbour will be the last curve
                        if (j == 0)
                            Denominator1.Add(d.Count - 1);
                        else
                            Denominator1.Add(j - 1);
                        double ProVal = ProductFunction(d, Denominator1);
                        denominator += ProVal;
                    }
                    Kapa = Numerator / denominator;
                    Mus.Add(Kapa);
                }
            }
            if (method == Blending_Method.Special_Side_Blending)
            {
                for (int i = 0; i < d.Count; i++)
                {
                    double Nu = 0;
                    double Numerator = 0;
                    double denominator = 0;
                    List<int> Numerator1 = new List<int>();
                    Numerator1.Add(i);
                    Numerator = ProductFunction(d, Numerator1);
                    for (int j = 0; j < d.Count; j++)
                    {
                        List<int> Denominator1 = new List<int>();
                        Denominator1.Add(j);
                        double ProVal = ProductFunction(d, Denominator1);
                        denominator += ProVal;
                    }
                    Nu = Numerator / denominator;
                    if (double.IsNaN(Nu))
                        Nu = 0;
                    Mus.Add(Nu);
                }
                //limit constraint for corner points in order to avoid singularity
                for (int i = 0; i < Mus.Count; i++)
                {
                    var side = i;
                    var Preside = i - 1;
                    //If boundary curve ID equals to 0,its previous neighbour will be the last curve
                    if (i == 0)
                        Preside = d.Count - 1;
                    if (d[side] < 0.0001 && d[Preside] < 0.0001)
                    {
                        Mus[side] = 1;
                        Mus[Preside] = 0;
                    }
                }
            }
            if (method == Blending_Method.Corner_Blending)
            {
                for (int i = 0; i < d.Count; i++)
                {
                    double Kapa = 0;
                    double Numerator = 0;
                    double denominator = 0;
                    List<int> Numerator1 = new List<int>();
                    List<int> Numerator2 = new List<int>();
                    Numerator1.Add(i);
                    //If boundary curve ID equals to 0,its previous neighbour will be the last curve
                    if (i == 0)
                        Numerator1.Add(d.Count - 1);
                    else
                        Numerator1.Add(i - 1);
                    for (int j = 0; j < d.Count; j++)
                    {
                        List<int> Denominator1 = new List<int>();
                        Denominator1.Add(j);
                        //If boundary curve ID equals to 0,its previous neighbour will be the last curve
                        if (j == 0)
                            Denominator1.Add(d.Count - 1);
                        else
                            Denominator1.Add(j - 1);
                        double ProVal = ProductFunction(d, Denominator1);
                        denominator += ProVal;
                    }
                    Kapa = Numerator / denominator;
                    Mus.Add(Kapa);
                }
            }
            return Mus;
        }
        private static double ProductFunction(List<double> d, List<int> n)
        {
            double Result = 0;
            List<double> SqreList = new List<double>();
            for (int i = 0; i < d.Count; i++)
            {
                bool eliminated = false;
                for (int j = 0; j < n.Count; j++)
                {
                    if (i == n[j])
                    {
                        eliminated = true;
                        break;
                    }
                }
                if (eliminated)
                    continue;
                SqreList.Add(Math.Pow(d[i], 2));
            }
            double Product = 1.0;
            for (int i = 0; i < SqreList.Count; i++)
            {
                Product *= SqreList[i];
            }
            Result = Product;
            return Result;
        }
    }
}






