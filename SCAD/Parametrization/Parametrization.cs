using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using Rhino.Input.Custom;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static SCAD.Extensions;
using static SCAD.HarmonicClass;

namespace SCAD
{
    public enum Parametrization_Method
    {
        //
        // Summary:
        //      Uses Side Blending functions
        IsoLine = 0,
        //
        // Summary:
        //      Uses Corner Blending functions
        MVC,
        //
        // Summary:
        //      Uses Special Side Blending functions
        Harmonic_C,
        //
        // Summary:
        //      Uses Radial distance functions were suggested by Charrot and Gregory, and these also work for non-regular domains. Input domains should be lines.
        Harmonic_Pt,
        RadialDistanceFunction
    }

    public class Parametrization
    {
        public Parametrization_Method method = 0;
        public List<Curve> domainC = new List<Curve>();
        public List<List<Point3d>> domainPt = new List<List<Point3d>>();
        public List<Point3d> domainL = new List<Point3d>();
        public List<HarmonicMap> HarmonicMapList_si = new List<HarmonicMap>();
        public List<HarmonicMap> HarmonicMapList_di = new List<HarmonicMap>();
        private Domain domain;


        public Parametrization() { }

        public Parametrization(Parametrization_Method method, List<Curve> domaincurves, Domain domain)
        {
            this.method = method;
            this.domainC = domaincurves;
            this.domain = domain;
            //HarmonicMapCreate_curve(domainC);
        }
        public Parametrization(Parametrization_Method method, Domain domain)
        {
            if (method == Parametrization_Method.Harmonic_Pt)
            {
                var newcurves = new List<List<Point3d>>();
                for (int i = 0; i < domain.Curves.Count; i++)
                {
                    NurbsCurve nb = domain.Curves[i].ToNurbsCurve();
                    List<Point3d> pts = new List<Point3d>();
                    for (int j = 0; j < nb.Points.Count; j++)
                        pts.Add(new Point3d(nb.Points[j].X, nb.Points[j].Y, 0));
                    newcurves.Add(pts);
                }

                this.method = method;
                this.domainPt = newcurves;
                this.domain = domain;
                HarmonicMapCreate_curve(domainPt);
            }
            else if (method == Parametrization_Method.RadialDistanceFunction)
            {
                this.method = method;
                this.domainL = domain.Vertices3d;
                this.domain = domain;
            }
        }

        public (List<double> si, List<double> di) GetPoint(double u, double v)
        {
            List<double> si = new List<double>();
            List<double> di = new List<double>();

            if (method == Parametrization_Method.Harmonic_C)
            {
                (si, di) = Harmonic(u, v, domainC.Count);
            }
            else if (method == Parametrization_Method.Harmonic_Pt)
            {
                (si, di) = Harmonic(u, v, domainPt.Count);
            }
            else if (method == Parametrization_Method.RadialDistanceFunction)
            {
                (si, di) = RadialDistanceFunction(u, v, domainL);
            }
            return (si: si, di: di);
        }

        private (List<double> si, List<double> di) Harmonic(double u, double v, int n)
        {
            List<double> si = new List<double>();
            List<double> di = new List<double>();

            int j;
            Point3d point = new Point3d(u, v, 0);
            for (int i = 0; i < n; i++)
            {
                ///
                bool success_si = Harmonic_eval(HarmonicMapList_si[i], point, out double result_si);//harmonic calculation
                ///
                bool success_di = Harmonic_eval(HarmonicMapList_di[i], point, out double result_di);//harmonic calculation
                ///

                //if (success_si == false) AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Harmonic function si!!");
                //if (success_di == false) AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Harmonic function di!!");
                if (success_si == false) System.Diagnostics.Debug.WriteLine("Error Harmonic function si!!");
                if (success_di == false) System.Diagnostics.Debug.WriteLine("Error Harmonic function di!!");

                si.Add(result_si);
                di.Add(result_di);
            }
            return (si: si, di: di);

        }

        private void HarmonicMapCreate_curve(List<List<Point3d>> curvepointlist)
        {
            HarmonicMapList_si = new List<HarmonicMap>();
            HarmonicMapList_di = new List<HarmonicMap>();
            List<(double, double)> output = new List<(double, double)>();
            int i_before, i_after, selected_curve = 2;
            bool biharmonic = true;
            double epsilon = 1.0e-5;
            int levels = 9;
            double value = 1;
            (Interval Ix, Interval Iy) = domain.Bounds;
            double[] minn = { Ix.Min, Iy.Min }, maxx = { Ix.Max, Iy.Max }; // domain polygon size

            //// Domain curve creation without "value" or "z"
            List<Point3d[]> DomainPolygonHarmonic_main = new List<Point3d[]> { };
            for (int i = 0; i < curvepointlist.Count; i++)
                DomainPolygonHarmonic_main.Add(curvepointlist[i].ToArray());

            for (int i = 0; i < DomainPolygonHarmonic_main.Count; i++)
            {
                i_before = IndexWrapper((i - 1), DomainPolygonHarmonic_main.Count);
                int i_beforebefore = IndexWrapper((i - 2), DomainPolygonHarmonic_main.Count);
                i_after = IndexWrapper((i + 1), DomainPolygonHarmonic_main.Count);

                ///***********
                ///harmonic map created for si -----------------------------------
                ///***********
                double[] min = minn, max = maxx; // domain polygon size
                HarmonicMap map_si;
                map_si = Harmonic_create(min, max, levels);

                ///Assigned Value on domain points, si
                List<Point3d[]> DomainPolygonHarmonic = new List<Point3d[]>();
                for (int ii = 0; ii < DomainPolygonHarmonic_main.Count; ii++)
                {
                    DomainPolygonHarmonic.Add(new Point3d[DomainPolygonHarmonic_main[ii].Length]);
                    if (ii == i_after)
                    {
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, value);
                    }
                    else if (ii == i_before)
                    {
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, 0);
                    }
                    else if (ii == i)
                    {
                        double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, stepsize * pt);
                    }
                    else if (ii == i_beforebefore)
                    {
                        double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y,
                                stepsize * ((DomainPolygonHarmonic[ii].Length - 1) - pt));
                    }
                    else
                    {
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, value);
                    }
                }
                ////

                for (int ii = 0; ii < DomainPolygonHarmonic.Count; ++ii)
                {
                    ///
                    ///middle point repeat, bunu hocanın örnek modeli yüzünden eklemiştim, 999 ise iptal etmişizdir.
                    if (ii == selected_curve)
                    {
                        List<Point3d> idle = DomainPolygonHarmonic[ii].ToList();
                        for (int xx = 0; xx < idle.Count; xx++)
                        {
                            if (xx == 3)
                            {
                                Point3d asd = new Point3d(idle[xx]);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                break;
                            }
                        }
                        DomainPolygonHarmonic[ii] = idle.ToArray();
                    }
                    ////
                    ///
                    Harmonic_add_curve(map_si, DomainPolygonHarmonic[ii], DomainPolygonHarmonic[ii].Length);
                }
                Harmonic_solve(map_si, epsilon, biharmonic);
                ///****
                HarmonicMapList_si.Add(map_si);
                ////
                ////--------------------------------------------------------------

                ///***********
                ///harmonic map created for di -----------------------------------
                ///***********
                //double[] min_di = { Ix.Min, Iy.Min }, max_di = { Ix.Max, Iy.Max };  // domain polygon size
                double[] min_di = minn, max_di = maxx; // domain polygon size
                HarmonicMap map_di;
                map_di = Harmonic_create(min_di, max_di, levels);

                ///Assigned Value on domain points, di
                DomainPolygonHarmonic = new List<Point3d[]>();
                for (int ii = 0; ii < DomainPolygonHarmonic_main.Count; ii++)
                {
                    DomainPolygonHarmonic.Add(new Point3d[DomainPolygonHarmonic_main[ii].Length]);
                    if (ii == i_after)
                    {
                        double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, stepsize * pt);
                    }
                    else if (ii == i_before)
                    {
                        double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, stepsize * ((DomainPolygonHarmonic[ii].Length - 1) - pt));
                    }
                    else if (ii == i)
                    {
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, 0);
                    }
                    else
                    {
                        for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
                            DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, value);
                    }
                }
                ////

                for (int ii = 0; ii < DomainPolygonHarmonic.Count; ++ii)
                {
                    ///
                    ///middle point repeat
                    if (ii == selected_curve)
                    {
                        List<Point3d> idle = DomainPolygonHarmonic[ii].ToList();
                        for (int xx = 0; xx < idle.Count; xx++)
                        {
                            if (xx == 3)
                            {
                                Point3d asd = new Point3d(idle[xx]);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                idle.Insert(xx, asd);
                                break;
                            }
                        }
                        DomainPolygonHarmonic[ii] = idle.ToArray();
                    }
                    ////
                    ///
                    Harmonic_add_curve(map_di, DomainPolygonHarmonic[ii], DomainPolygonHarmonic[ii].Length);
                }
                Harmonic_solve(map_di, epsilon, biharmonic);
                ///*****
                HarmonicMapList_di.Add(map_di);
                ////
                ////-------------------------------------------------------------
            }

        }

        //private void HarmonicMapCreate_curve(List<Curve> domaincurves)
        //{
        //    HarmonicMapList_si = new List<HarmonicMap>();
        //    HarmonicMapList_di = new List<HarmonicMap>();
        //    List<(double, double)> output = new List<(double, double)>();
        //    int i_before, i_after, selected_curve = 999;
        //    bool biharmonic = false;
        //    double epsilon = 1.0e-5;
        //    int levels = 3;
        //    double value = 1;
        //    (Interval Ix, Interval Iy) = domain.Bounds;
        //    double[] minn = { Ix.Min, Iy.Min }, maxx = { Ix.Max, Iy.Max }; // domain polygon size

        //    ///// Extraction Rhino curve control points for harmonic function
        //    List<List<Point3d>> curvepointlist = new List<List<Point3d>>();
        //    for (int i = 0; i < domaincurves.Count; i++)
        //    {
        //        List<Point3d> pts = new List<Point3d>();
        //        foreach (var item in domaincurves[i].ToNurbsCurve().Points)
        //        {
        //            pts.Add(item.Location);
        //        }
        //        curvepointlist.Add(pts);
        //    }


        //    //// Domain curve creation without "value" or "z"
        //    List<Point3d[]> DomainPolygonHarmonic_main = new List<Point3d[]> { };
        //    for (int i = 0; i < curvepointlist.Count; i++)
        //    {
        //        DomainPolygonHarmonic_main.Add(curvepointlist[i].ToArray());
        //    }

        //    for (int i = 0; i < DomainPolygonHarmonic_main.Count; i++)
        //    {
        //        i_before = IndexWrapper((i - 1), DomainPolygonHarmonic_main.Count);
        //        int i_beforebefore = IndexWrapper((i - 2), DomainPolygonHarmonic_main.Count);
        //        i_after = IndexWrapper((i + 1), DomainPolygonHarmonic_main.Count);

        //        ///***********
        //        ///harmonic map created for si -----------------------------------
        //        ///***********
        //        double[] min = minn, max = maxx; // domain polygon size
        //        HarmonicMap map_si;
        //        map_si = harmonic_create(min, max, levels);

        //        ///Assigned Value on domain points, si
        //        List<Point3d[]> DomainPolygonHarmonic = new List<Point3d[]>();
        //        for (int ii = 0; ii < DomainPolygonHarmonic_main.Count; ii++)
        //        {
        //            DomainPolygonHarmonic.Add(new Point3d[DomainPolygonHarmonic_main[ii].Length]);
        //            if (ii == i_after)
        //            {
        //                for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
        //                {
        //                    DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, value);
        //                }
        //            }
        //            else if (ii == i_before)
        //            {
        //                for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
        //                {
        //                    DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, 0);
        //                }
        //            }
        //            else if (ii == i)
        //            {
        //                double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
        //                for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
        //                {
        //                    DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, stepsize * pt);
        //                }
        //            }
        //            else if (ii == i_beforebefore)
        //            {
        //                double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
        //                for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
        //                {
        //                    DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, stepsize * ((DomainPolygonHarmonic[ii].Length - 1) - pt));
        //                }
        //            }
        //            else
        //            {
        //                for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
        //                {
        //                    DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, value);
        //                }
        //            }
        //        }
        //        ////

        //        for (int ii = 0; ii < DomainPolygonHarmonic.Count; ++ii)
        //        {
        //            ///
        //            ///middle point repeat
        //            if (ii == selected_curve)
        //            {
        //                List<Point3d> idle = DomainPolygonHarmonic[ii].ToList();
        //                for (int xx = 0; xx < idle.Count; xx++)
        //                {
        //                    if (xx == 3)
        //                    {
        //                        Point3d asd = new Point3d(idle[xx]);
        //                        idle.Insert(xx, asd);
        //                        idle.Insert(xx, asd);
        //                        idle.Insert(xx, asd);
        //                        idle.Insert(xx, asd);
        //                        //idle.Insert(xx, asd);
        //                        break;
        //                    }
        //                }
        //                DomainPolygonHarmonic[ii] = idle.ToArray();
        //            }
        //            ////
        //            ///
        //            harmonic_add_curve(map_si, DomainPolygonHarmonic[ii], DomainPolygonHarmonic[ii].Length);
        //        }
        //        harmonic_solve(map_si, 1.0e-5, false);
        //        ///****
        //        HarmonicMapList_si.Add(map_si);
        //        ////
        //        ////--------------------------------------------------------------

        //        ///***********
        //        ///harmonic map created for di -----------------------------------
        //        ///***********
        //        double[] min_di = minn, max_di = maxx; // domain polygon size
        //        HarmonicMap map_di;
        //        map_di = harmonic_create(min_di, max_di, levels);

        //        ///Assigned Value on domain points, di
        //        DomainPolygonHarmonic = new List<Point3d[]>();
        //        for (int ii = 0; ii < DomainPolygonHarmonic_main.Count; ii++)
        //        {
        //            DomainPolygonHarmonic.Add(new Point3d[DomainPolygonHarmonic_main[ii].Length]);
        //            if (ii == i_after)
        //            {
        //                double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
        //                for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
        //                {
        //                    DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, stepsize * pt);
        //                }
        //            }
        //            else if (ii == i_before)
        //            {
        //                double stepsize = value / (DomainPolygonHarmonic[ii].Length - 1);
        //                for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
        //                {
        //                    DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, stepsize * ((DomainPolygonHarmonic[ii].Length - 1) - pt));
        //                }
        //            }
        //            else if (ii == i)
        //            {
        //                for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
        //                {
        //                    DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, 0);
        //                }
        //            }
        //            else
        //            {
        //                for (int pt = 0; pt < DomainPolygonHarmonic[ii].Length; pt++)
        //                {
        //                    DomainPolygonHarmonic[ii][pt] = new Point3d(DomainPolygonHarmonic_main[ii][pt].X, DomainPolygonHarmonic_main[ii][pt].Y, value);
        //                }
        //            }
        //        }
        //        ////

        //        for (int ii = 0; ii < DomainPolygonHarmonic.Count; ++ii)
        //        {
        //            ///
        //            ///middle point repeat
        //            if (ii == selected_curve)
        //            {
        //                List<Point3d> idle = DomainPolygonHarmonic[ii].ToList();
        //                for (int xx = 0; xx < idle.Count; xx++)
        //                {
        //                    if (xx == 3)
        //                    {
        //                        Point3d asd = new Point3d(idle[xx]);
        //                        idle.Insert(xx, asd);
        //                        idle.Insert(xx, asd);
        //                        idle.Insert(xx, asd);
        //                        idle.Insert(xx, asd);
        //                        idle.Insert(xx, asd);
        //                        break;
        //                    }
        //                }
        //                DomainPolygonHarmonic[ii] = idle.ToArray();
        //            }
        //            ////
        //            ///
        //            harmonic_add_curve(map_di, DomainPolygonHarmonic[ii], DomainPolygonHarmonic[ii].Length);
        //        }
        //        harmonic_solve(map_di, 1.0e-5, false);
        //        ///*****
        //        HarmonicMapList_di.Add(map_di);
        //        ////
        //        ////-------------------------------------------------------------
        //    }

        //}

        /// <summary>
        /// Radial distance functions were suggested by Charrot and Gregory, and these also work for non-regular domains.
        /// </summary>
        /// <param name="u"> u value</param>
        /// <param name="v"> v value</param>
        /// <param name="domainlines"> domain lines</param>
        /// <param name="si_di"> local parametrization distance from each curve, first item is s value, second item is d value</param>
        private (List<double> si, List<double> di) RadialDistanceFunction(double u, double v, List<Point3d> domainlines)
        {
            List<double> si = new List<double>();
            List<double> di = new List<double>();
            int i_before, i_after, i_afterafter;

            int n = domainlines.Count;
            List<Point3d> domainpoints = new List<Point3d>();
            //foreach (var item in domainlines)
            //{
            //    domainpoints.Add(item.From);
            //}
            domainpoints = domainlines;
            for (int i = 0; i < n; i++)
            {

                i_before = IndexWrapper((i - 1), n);
                i_after = IndexWrapper((i + 1), n);
                i_afterafter = IndexWrapper((i + 2), n);

                Line i_beforeline = new Line(domainpoints[i], domainpoints[i_before]);
                Line i_afterline = new Line(domainpoints[i_after], domainpoints[i_afterafter]);

                bool A = Intersection.LineLine(i_beforeline, i_afterline, out double a, out double b, 0.001, false);
                Point3d c_i = i_beforeline.PointAt(a);
                Point3d c_i_0 = i_afterline.PointAt(b);


                Line null_line = new Line(new Point3d(u, v, 0), c_i);
                Line null_line2 = new Line(domainpoints[i], domainpoints[i_after]);
                bool B = Intersection.LineLine(null_line, null_line2, out double a1, out double b2, 0.001, true);
                Point3d e_i = null_line.PointAt(a1);
                Point3d e_i0 = null_line2.PointAt(b2);

                if (!A || !B || (c_i - c_i_0).Length > 0.000001 || (e_i - e_i0).Length > 0.000001)
                    System.Diagnostics.Debug.WriteLine("Radial distance function domain line Intersection problem!!");
                //AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Radial distance function domain line Intersection problem!!");

                si.Add((e_i - domainpoints[i]).Length / (domainpoints[i_after] - domainpoints[i]).Length);
                di.Add((new Point3d(u, v, 0) - e_i).Length);
            }
            return (si: si, di: di);
        }
    }

    public static class HarmonicClass
    {
        static int MIN_LEVEL = 3;
        static double MAX(double a, double b) => (a > b) ? a : b;

        public struct GridValue
        {
            public bool Boundary { get; set; }
            public double Value { get; set; }
        };

        public struct HarmonicMap
        {
            public int Levels { get; set; }
            public int Size { get; set; }
            public GridValue[] Grid { get; set; }
            public double[] Offset { get; set; }
            public double Scaling { get; set; }
        };


        public static HarmonicMap Harmonic_create(double[] min, double[] max, int levels)
        {
            /* Create the grid */
            int n = (int)Math.Pow(2, levels);
            HarmonicMap map = new HarmonicMap();
            map.Offset = new double[2];
            map.Levels = levels;
            map.Size = n;
            map.Grid = new GridValue[n * n];

            /* Add a margin of 2.5% on all sides */
            double length = MAX(max[0] - min[0], max[1] - min[1]);
            map.Offset[0] = min[0] - length * 0.025;
            map.Offset[1] = min[1] - length * 0.025;
            map.Scaling = (double)n / length / 1.05;

            /* Initialize cells */
            for (int i = 0; i < n * n; ++i)
            {
                map.Grid[i].Boundary = false;
                map.Grid[i].Value = 0.0;
            }
            return map;
        }

        public static void SolveHarmonic(GridValue[] grid, int n, double epsilon)
        {
            double change;
            do
            {
                change = 0.0;
                int count = 0, index = n + 1;
                for (int j = 1, n_1 = n - 1; j < n_1; ++j)
                {
                    for (int i = 1; i < n_1; ++i, ++index)
                        if (!grid[index].Boundary)
                        {
                            double value = 0.0;
                            value += grid[index - n].Value;
                            value += grid[index - 1].Value;
                            value += grid[index + n].Value;
                            value += grid[index + 1].Value;
                            value /= 4.0;
                            change += Math.Abs(grid[index].Value - value);
                            grid[index].Value = value;
                            ++count;
                        }
                    index += 2;
                }
                /* Boundary cases: not handled [the result is the same] */
                change /= (double)count;
            } while (change > epsilon);
        }

        public static void SolveBiharmonic(GridValue[] grid, int n, double epsilon)
        {
            double change;
            do
            {
                change = 0.0;
                int count = 0, n_2 = n - 2, n2 = n * 2, index = n2 + 2;
                for (int j = 2; j < n_2; ++j)
                {
                    for (int i = 2; i < n_2; ++i, ++index)
                        if (!grid[index].Boundary)
                        {
                            double value = 0.0;
                            value += grid[index - n].Value;
                            value += grid[index - 1].Value;
                            value += grid[index + n].Value;
                            value += grid[index + 1].Value;
                            value *= 4.0;
                            value -= grid[index - n - 1].Value;
                            value -= grid[index - n + 1].Value;
                            value -= grid[index + n - 1].Value;
                            value -= grid[index + n + 1].Value;
                            value *= 2.0;
                            value -= grid[index - n2].Value;
                            value -= grid[index - 2].Value;
                            value -= grid[index + n2].Value;
                            value -= grid[index + 2].Value;
                            value /= 20.0;
                            change += Math.Abs(grid[index].Value - value);
                            grid[index].Value = value;
                            ++count;
                        }
                    index += 4;
                }
                /* Boundary cases */
                for (int j = 0; j < n; ++j)
                    for (int i = 0; i < n; ++i)
                    {
                        if (j >= 2 && j < n_2 && i >= 2 && i < n_2)
                            continue;
                        index = j * n + i;
                        if (grid[index].Boundary)
                            continue;
                        double weight = 0.0;
                        double value = 0.0;
                        if (i > 0)
                        {
                            int k = index - 1;
                            weight += 2.0;
                            value += grid[k].Value * 2.0;
                            int neighbors = 4;
                            if (i == 1)
                                --neighbors;
                            if (j == 0)
                                --neighbors;
                            if (j == n - 1)
                                --neighbors;
                            double w = 1.0 / (double)neighbors;
                            if (i > 1)
                            {
                                weight -= w;
                                value -= grid[k - 1].Value * w;
                            }
                            if (j > 0)
                            {
                                weight -= w;
                                value -= grid[k - n].Value * w;
                            }
                            if (j < n - 1)
                            {
                                weight -= w;
                                value -= grid[k + n].Value * w;
                            }
                        }
                        if (i < n - 1)
                        {
                            int k = index + 1;
                            weight += 2.0;
                            value += grid[k].Value * 2.0;
                            int neighbors = 4;
                            if (i == n - 2)
                                --neighbors;
                            if (j == 0)
                                --neighbors;
                            if (j == n - 1)
                                --neighbors;
                            double w = 1.0 / (double)neighbors;
                            if (i < n - 2)
                            {
                                weight -= w;
                                value -= grid[k + 1].Value * w;
                            }
                            if (j > 0)
                            {
                                weight -= w;
                                value -= grid[k - n].Value * w;
                            }
                            if (j < n - 1)
                            {
                                weight -= w;
                                value -= grid[k + n].Value * w;
                            }
                        }
                        if (j > 0)
                        {
                            int k = index - n;
                            weight += 2.0;
                            value += grid[k].Value * 2.0;
                            int neighbors = 4;
                            if (i == 0)
                                --neighbors;
                            if (i == n - 1)
                                --neighbors;
                            if (j == 1)
                                --neighbors;
                            double w = 1.0 / (double)neighbors;
                            if (i > 0)
                            {
                                weight -= w;
                                value -= grid[k - 1].Value * w;
                            }
                            if (i < n - 1)
                            {
                                weight -= w;
                                value -= grid[k + 1].Value * w;
                            }
                            if (j > 1)
                            {
                                weight -= w;
                                value -= grid[k - n].Value * w;
                            }
                        }
                        if (j < n - 1)
                        {
                            int k = index + n;
                            weight += 2.0;
                            value += grid[k].Value * 2.0;
                            int neighbors = 4;
                            if (i == 0)
                                --neighbors;
                            if (i == n - 1)
                                --neighbors;
                            if (j == n - 2)
                                --neighbors;
                            double w = 1.0 / (double)neighbors;
                            if (i > 0)
                            {
                                weight -= w;
                                value -= grid[k - 1].Value * w;
                            }
                            if (i < n - 1)
                            {
                                weight -= w;
                                value -= grid[k + 1].Value * w;
                            }
                            if (j < n - 2)
                            {
                                weight -= w;
                                value -= grid[k + n].Value * w;
                            }
                        }
                        value /= weight;
                        change += Math.Abs(grid[index].Value - value);
                        grid[index].Value = value;
                        ++count;
                    }
                change /= (double)count;
            } while (change > epsilon);
        }

        public static void Solve(GridValue[] grid, int level, double epsilon, bool biharmonic)
        {
            int n = (int)Math.Pow(2, level);
            if (level > MIN_LEVEL)
            {
                /* Generate a coarser grid and solve that first to get good starting values */
                int level1 = level - 1, n1 = (int)Math.Pow(2, level1);
                GridValue[] grid1 = new GridValue[n1 * n1];
                for (int i = 0; i < n1; ++i)
                    for (int j = 0; j < n1; ++j)
                    {
                        grid1[j * n1 + i].Value = 0.0;
                        int count = 0;
                        if (grid[2 * j * n + 2 * i].Boundary)
                        {
                            ++count;
                            grid1[j * n1 + i].Value += grid[2 * j * n + 2 * i].Value;
                        }
                        if (grid[2 * j * n + 2 * i + 1].Boundary)
                        {
                            ++count;
                            grid1[j * n1 + i].Value += grid[2 * j * n + 2 * i + 1].Value;
                        }
                        if (grid[(2 * j + 1) * n + 2 * i].Boundary)
                        {
                            ++count;
                            grid1[j * n1 + i].Value += grid[(2 * j + 1) * n + 2 * i].Value;
                        }
                        if (grid[(2 * j + 1) * n + 2 * i + 1].Boundary)
                        {
                            ++count;
                            grid1[j * n1 + i].Value += grid[(2 * j + 1) * n + 2 * i + 1].Value;
                        }
                        if (count > 0)
                        {
                            grid1[j * n1 + i].Boundary = true;
                            grid1[j * n1 + i].Value /= (double)count;
                        }
                        else
                            grid1[j * n1 + i].Boundary = false;
                    }
                Solve(grid1, level1, epsilon, biharmonic);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        if (!grid[j * n + i].Boundary)
                            grid[j * n + i].Value = grid1[(j / 2) * n1 + i / 2].Value;
                grid1 = null;
            }

            /* Solve by iteration */
            if (biharmonic)
                SolveBiharmonic(grid, n, epsilon);
            else
                SolveHarmonic(grid, n, epsilon);
        }

        public static void Harmonic_add_point(HarmonicMap map, Point3d point)
        {
            int n = map.Size;
            int x = (int)Math.Round((point[0] - map.Offset[0]) * map.Scaling);
            int y = (int)Math.Round((point[1] - map.Offset[1]) * map.Scaling);
            map.Grid[y * n + x].Boundary = true;
            map.Grid[y * n + x].Value = point[2];
        }

        public static void Harmonic_add_line(HarmonicMap map, Point3d from, Point3d to)
        {
            int n = map.Size;
            int x0 = (int)Math.Round((from[0] - map.Offset[0]) * map.Scaling);
            int y0 = (int)Math.Round((from[1] - map.Offset[1]) * map.Scaling);
            double v0 = from[2];
            int x1 = (int)Math.Round((to[0] - map.Offset[0]) * map.Scaling);
            int y1 = (int)Math.Round((to[1] - map.Offset[1]) * map.Scaling);
            double v1 = to[2];
            int dx = Math.Abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
            int dy = Math.Abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
            int err = (dx > dy ? dx : -dy) / 2, e2;
            if (err == 0)
            {
                map.Grid[y0 * n + x0].Boundary = true;
                map.Grid[y0 * n + x0].Value = v0;
                map.Grid[y1 * n + x1].Boundary = true;
                map.Grid[y1 * n + x1].Value = v1;
                return;
            }
            while (true)
            {
                double ratio;             /* linear interpolation along the sides */
                if (err > 0)
                    ratio = (double)Math.Abs(x1 - x0) / (double)dx;
                else
                    ratio = (double)Math.Abs(y1 - y0) / (double)dy;
                map.Grid[y0 * n + x0].Boundary = true;
                map.Grid[y0 * n + x0].Value = v0 * ratio + v1 * (1.0 - ratio);
                if (x0 == x1 && y0 == y1) break;
                e2 = err;
                if (e2 > -dx) { err -= dy; x0 += sx; }
                if (e2 < dy) { err += dx; y0 += sy; }
            }
        }

        public static void Harmonic_add_curve(HarmonicMap map, Point3d[] points, int n)
        {
            double tmp;
            int resolution = map.Size;
            double[] coeff = new double[n];
            Point3d from, to;
            from = new Point3d(points[0].X, points[0].Y, points[0].Z);
            for (int i = 1; i <= resolution; ++i)
            {
                double u = (double)i / resolution;
                /* Compute Bernstein polynomials */
                coeff[0] = 1.0;
                for (int j = 1; j < n; ++j)
                {
                    double saved = 0.0;
                    for (int k = 0; k < j; ++k)
                    {
                        tmp = coeff[k];
                        coeff[k] = saved + tmp * (1.0 - u);
                        saved = tmp * u;
                    }
                    coeff[j] = saved;
                }
                /* Evaluate the curve */
                to = new Point3d(0.0, 0.0, 0.0);
                for (int j = 0; j < n; ++j)
                {
                    to += points[j] * coeff[j];
                }
                /* Draw a segment */
                Harmonic_add_line(map, from, to);
                /* Swap from & to */
                (from, _) = (to, from);
            }
            coeff = null;
        }

        public static void Harmonic_solve(HarmonicMap map, double epsilon, bool biharmonic)
        {
            Solve(map.Grid, map.Levels, epsilon, biharmonic);
        }

        public static bool Inside_map(HarmonicMap map, int i, int j)
        {
            return i >= 0 && j >= 0 && i < map.Size && j < map.Size;
        }

        public static bool Harmonic_eval(HarmonicMap map, Point3d point, out double value)
        {
            int n = map.Size;
            double x = (point[0] - map.Offset[0]) * map.Scaling;
            double y = (point[1] - map.Offset[1]) * map.Scaling;
            int i = (int)Math.Round(x), j = (int)Math.Round(y);

            if (!(Inside_map(map, i, j) &&
                  Inside_map(map, i, j + 1) &&
                  Inside_map(map, i + 1, j) &&
                  Inside_map(map, i + 1, j + 1)))
            {
                value = 0;
                return false;               /* The point is outside the region */
            }

            value = map.Grid[j * n + i].Value * (1.0 - y + j) * (1.0 - x + i);
            value += map.Grid[(j + 1) * n + i].Value * (y - j) * (1.0 - x + i);
            value += map.Grid[j * n + i + 1].Value * (1.0 - y + j) * (x - i);
            value += map.Grid[(j + 1) * n + i + 1].Value * (y - j) * (x - i);

            return true;
        }

        public static void Harmonic_write_ppm(HarmonicMap map, string filename)
        {
            int n = map.Size;
            StreamWriter f = new StreamWriter(filename);

            f.Write("P3\n%zu %zu\n255\n", n, n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                    if (map.Grid[j * n + i].Boundary)
                        f.Write("255 0 0 ");
                    else
                        f.Write("0 0 " + ((int)Math.Round(map.Grid[j * n + i].Value * 255.0)).ToString("N6"));
                //f.Write("0 0 %d ", (int)Math.Round(map.grid[j * n + i].value * 255.0));
                f.Write("\n");
            }
            f.Close();
        }

        public static void Harmonic_free(HarmonicMap map)
        {
            map.Grid = null;
            map = new HarmonicMap();
        }

    }

    public class MVCClass
    {
        private List<double> w_i = new List<double>();
        private List<double> s_i = new List<double>();

        public List<(double, double)> MVC(List<Point3d> domainPolygon, double u, double v)
        {
            var pt = new Point3d(u, v, 0);
            MVCweights(domainPolygon, pt, out w_i);

            var d_i = w_i.Select(x => 1 - Math.Abs(x)).ToList();
            for (int i = 0; i < w_i.Count; i++)
            {
                if (w_i[i] == 0 && w_i[(i + 1) % w_i.Count] == 0)
                    s_i.Add(1);
                else
                    s_i.Add((double)d_i[i] / (double)(d_i[i] + d_i[(i + 1) % d_i.Count]));
            }
            for (int i = 0; i < d_i.Count; i++)
            {
                d_i[i] = (1 - Math.Abs(w_i[i] + (w_i[(i + 1) % w_i.Count] - w_i[i])) * s_i[i]);
            }
            List<(double, double)> Distance = new List<(double, double)>();

            for (int i = 0; i < d_i.Count; i++)
                Distance.Add((s_i[i], d_i[i]));

            return Distance;
        }

        private static int MVCweights(List<Point3d> cageCoords, Point3d queryCoord, out List<double> baryCoords)
        {
            int nSize = cageCoords.Count;
            Debug.Assert(nSize != 0);

            double dx, dy;

            List<Vector3d> s = new List<Vector3d>(nSize);
            baryCoords = new List<double>(nSize);

            for (int i = 0; i < nSize; i++)
            {
                dx = cageCoords[i].X - queryCoord.X;
                dy = cageCoords[i].Y - queryCoord.Y;
                s.Add(new Vector3d(dx, dy, 0));
            }

            for (int i = 0; i < nSize; i++)
                baryCoords.Add(0.0);

            int ip, im;      // (i+1) and (i-1)
            double ri, rp, Ai, Di, dl, mu;  // Distance
            double eps = 10.0 * 2.22507e-308;

            // First check if any coordinates close to the cage point or
            // lie on the cage boundary. These are special cases.
            for (int i = 0; i < nSize; i++)
            {
                ip = (i + 1) % nSize;
                ri = Math.Sqrt(s[i].X * s[i].X + s[i].Y * s[i].Y);
                Ai = 0.5 * (s[i].X * s[ip].Y - s[ip].X * s[i].Y);
                Di = s[ip].X * s[i].X + s[ip].Y * s[i].Y;
                if (ri <= eps)
                {
                    baryCoords[i] = 1.0;
                    return 0;
                }
                else if (Math.Abs(Ai) <= 0 && Di < 0.0)
                {
                    dx = cageCoords[ip].X - cageCoords[i].X;
                    dy = cageCoords[ip].Y - cageCoords[i].Y;
                    dl = Math.Sqrt(dx * dx + dy * dy);
                    Debug.Assert(dl > eps);
                    dx = queryCoord.X - cageCoords[i].X;
                    dy = queryCoord.Y - cageCoords[i].Y;
                    mu = Math.Sqrt(dx * dx + dy * dy) / dl;
                    Debug.Assert(mu >= 0.0 && mu <= 1.0);
                    baryCoords[i] = 1.0 - mu;
                    baryCoords[ip] = mu;
                    return 0;
                }
            }

            // Page #12, from the paper
            List<double> tanalpha = new List<double>(nSize); // tan(alpha/2)
            for (int i = 0; i < nSize; i++)
            {
                ip = (i + 1) % nSize;
                im = (nSize - 1 + i) % nSize;
                ri = Math.Sqrt(s[i].X * s[i].X + s[i].Y * s[i].Y);
                rp = Math.Sqrt(s[ip].X * s[ip].X + s[ip].Y * s[ip].Y);
                Ai = 0.5 * (s[i].X * s[ip].Y - s[ip].X * s[i].Y);
                Di = s[ip].X * s[i].X + s[ip].Y * s[i].Y;
                tanalpha.Add((ri * rp - Di) / (2.0 * Ai));
            }

            // Equation #11, from the paper
            double wi, wsum = 0.0;
            for (int i = 0; i < nSize; i++)
            {
                im = (nSize - 1 + i) % nSize;
                ri = Math.Sqrt(s[i].X * s[i].X + s[i].Y * s[i].Y);
                wi = 2.0 * (tanalpha[i] + tanalpha[im]) / ri;
                wsum += wi;
                baryCoords[i] = wi;
            }

            if (Math.Abs(wsum) > 0.0)
            {
                for (int i = 0; i < nSize; i++)
                    baryCoords[i] /= wsum;
            }

            return 0;

        }
    }

    public class IsoLineClass
    {
        private List<Curve> Isolines(List<Line> m_DPolygonLines, List<List<Curve>> IsolinesList)
        {
            IsolinesList.Clear();
            var Pnum = m_DPolygonLines.Count;
            for (int i = 0; i < Pnum; i++)
            {
                var lft = new Line();

                if (i == 0)
                    lft = m_DPolygonLines[Pnum - 1];
                else if (i == 2 || i == 5)
                    lft = m_DPolygonLines[(i + 1) % Pnum];
                else
                    lft = m_DPolygonLines[i - 1];


                var Left = lft.ToNurbsCurve();
                var Rght = new List<NurbsCurve>();
                if (i == 2 || i == 5)
                {
                    Rght.Add(m_DPolygonLines[(i - 1) % Pnum].ToNurbsCurve());
                    Rght.Add(m_DPolygonLines[(i - 2) % Pnum].ToNurbsCurve());
                    if (i == 2)
                        Rght.Add(m_DPolygonLines[Pnum - 1].ToNurbsCurve());
                    else
                        Rght.Add(m_DPolygonLines[(i - 3) % Pnum].ToNurbsCurve());
                }
                else
                {
                    Rght.Add(m_DPolygonLines[(i + 1) % Pnum].ToNurbsCurve());
                    Rght.Add(m_DPolygonLines[(i + 2) % Pnum].ToNurbsCurve());
                    Rght.Add(m_DPolygonLines[(i + 3) % Pnum].ToNurbsCurve());

                }


                var Right = Curve.JoinCurves(Rght);
                if (i == 0)
                {
                    //var rev2 = Left.Reverse();
                    var rev = Right[0].Reverse();
                }
                else if (i == 1)
                {
                    var rev2 = Left.Reverse();
                    //var rev = Right[0].Reverse();
                }
                else if (i == 2)
                {
                    var rev2 = Left.Reverse();
                    //var rev = Right[0].Reverse();
                }
                else if (i == 3)
                {
                    //var rev2 = Left.Reverse();
                    var rev = Right[0].Reverse();
                }
                else if (i == 4)
                {
                    var rev2 = Left.Reverse();
                    //var rev = Right[0].Reverse();
                }
                else if (i == 5)
                {
                    var rev2 = Left.Reverse();
                    //var rev = Right[0].Reverse();
                }

                IsolinesList.Add(Curve.CreateTweenCurves(Left, Right[0], 30, 0.0001).ToList());
                IsolinesList[i].Insert(0, Left);
                IsolinesList[i].Add(Right[0]);
                if (i == 2 || i == 5)
                    IsolinesList[i].Reverse();
                if (i == 0 || i == 2 || i == 3 || i == 5)
                    IsolinesList[i].ForEach(x => x.Reverse());
            }
            return IsolinesList[5];
        }
    }


}
