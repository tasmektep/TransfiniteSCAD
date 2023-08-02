using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;

namespace Grasshopper2
{
    public class Transfinite2 : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Transfinite2 class.
        /// </summary>
        public Transfinite2()
          : base("Transfinite2", "Transfinite2",
              "Transfinite2",
              "Category", "Transfinite2")
        {
        }

        private static Sphere m_Sph;
        private static List<Curve> m_Crvs;
        private static List<Line> m_RibTans;
        private static List<Line> m_Domain;
        private static Mesh m_DomainMesh;
        private static List<Brep> m_trimmedSphere;
        private static double m_RibScl;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddBrepParameter("Sphere", "Sphere", "Sphere", GH_ParamAccess.item);
            pManager.AddCurveParameter("Curves", "Curves", "Curves", GH_ParamAccess.list);
            pManager.AddLineParameter("RibTans", "RibTans", "RibTans", GH_ParamAccess.list);
            pManager.AddLineParameter("Domain", "Domain", "Domain", GH_ParamAccess.list);
            pManager.AddMeshParameter("Mesh", "Mesh", "Mesh", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            m_RibScl = 0.5;
            SetRibbons();
            SetDomain();

            SetMesh();
            //SetKatoSurface();
            //SetGregorySurface();
            SetMidPointPatch();

            DA.SetData(0, m_trimmedSphere[0]);
            DA.SetDataList(1, m_Crvs);
            DA.SetDataList(2, m_RibTans);
            DA.SetDataList(3, m_Domain);
            DA.SetData(4, m_DomainMesh);
        }

        private void SetMidPointPatch()
        {
            double u, v;
            List<Point3d> DomainPts = new List<Point3d>();
            List<double> s;
            List<double> d;
            Point3d Pt, I, Center2D;
            List<double> B;
            int k;
            Point3d P0;
            double B0;

            for (int i = 0; i < m_Domain.Count; i++)
                DomainPts.Add(m_Domain[i].From);

            Center2D = Point3d.Origin;
            for (int i = 0; i < m_Domain.Count; i++)
                Center2D = Center2D + m_Domain[i].From;
            Center2D = Center2D / m_Domain.Count;
            (s, d) = RadialDistanceFunction(Center2D.X, Center2D.Y, DomainPts);

            P0 = Point3d.Origin;
            for (int i = 0; i < m_Domain.Count; i++)
            {
                k = (i - 1 + s.Count) % s.Count;
                P0 = P0 + ComputeI(i, s.Count, s[i], s[k]);
            }
            P0 = P0 / m_Domain.Count;

            for (int i = 0; i < m_DomainMesh.Vertices.Count; i++)
            {
                u = m_DomainMesh.Vertices[i].X;
                v = m_DomainMesh.Vertices[i].Y;
                (s, d) = RadialDistanceFunction(u, v, DomainPts);
                B = GetMPBlending(s, d);

                B0 = 0;
                Pt = Point3d.Origin;
                for (int j = 0; j < s.Count; j++)
                {
                    k = (j - 1 + s.Count) % s.Count;
                    I = ComputeI(j, s.Count, s[j], s[k]);

                    Pt = Pt + I * B[j];
                    B0 = B0 + B[j];
                }
                B0 = 1 - B0;
                Pt = Pt + P0 * B0;

                m_DomainMesh.Vertices[i] = (Point3f)Pt;
            }
        }

        private Point3d ComputeI(int i, int sCount, double si, double si_1)
        {
            List<Point3d> DomainPts = new List<Point3d>();
            Point3d Ri_1, Ri, Q, I;
            double tIMin, tIMax, tIMax_1, tI, tI_1;
            Vector3d vecI, vecI_1, W;
            int i_1;

            tIMin = m_Crvs[i].Domain.Min;
            tIMax = m_Crvs[i].Domain.Max;
            tI = tIMin + si * (tIMax - tIMin);

            i_1 = (i - 1 + sCount) % sCount;
            tI_1 = m_Crvs[i_1].Domain.Min + si_1 * (m_Crvs[i_1].Domain.Max - m_Crvs[i_1].Domain.Min);

            tIMax_1 = m_Crvs[i_1].Domain.Max;
            vecI = GetRibbonTangent(i, tIMin);
            vecI_1 = GetRibbonTangent(i_1, tIMax_1);
            W = Vector3d.CrossProduct(vecI_1, vecI);
            Q = m_Crvs[i].PointAt(tIMin) + Gamma(1 - si_1) * vecI
                + Gamma(si) * vecI_1
                + Gamma(si) * Gamma(1 - si_1) * W;

            Ri_1 = m_Crvs[i_1].PointAt(tI_1) + si * GetRibbonTangent(i_1, tI_1);
            Ri = m_Crvs[i].PointAt(tI) + (1 - si_1) * GetRibbonTangent(i, tI);
            I = (Point3d)(Ri_1 + Ri - Q);

            return I;
        }

        private void SetGregorySurface()
        {
            double u, v;
            List<Point3d> DomainPts = new List<Point3d>();
            List<double> s;
            List<double> d;
            Point3d Pt, I;
            List<double> B;
            int k;

            for (int i = 0; i < m_Domain.Count; i++)
                DomainPts.Add(m_Domain[i].From);

            for (int i = 0; i < m_DomainMesh.Vertices.Count; i++)
            {
                u = m_DomainMesh.Vertices[i].X;
                v = m_DomainMesh.Vertices[i].Y;
                (s, d) = RadialDistanceFunction(u, v, DomainPts);
                B = GetBlending(d, 2);

                Pt = Point3d.Origin;
                for (int j = 0; j < s.Count; j++)
                {
                    k = (j - 1 + s.Count) % s.Count;
                    I = ComputeI(j, s.Count, s[j], s[k]);

                    Pt = Pt + I * B[j];
                }

                m_DomainMesh.Vertices[i] = (Point3f)Pt;
            }
        }

        private List<double> GetCornerBlending(List<double> d)
        {
            List<double> Ds = new List<double>();
            int j;
            double D;
            double Sum = 0;

            for (int i = 0; i < d.Count; i++)
            {
                j = (i - 1 + d.Count) % d.Count;
                Sum = Sum + 1.0 / (d[i] * d[i] * d[j] * d[j]);
            }

            for (int i = 0; i < d.Count; i++)
            {
                j = (i - 1 + d.Count) % d.Count;
                D = 1.0 / (d[i] * d[i] * d[j] * d[j] * Sum);

                Ds.Add(D);
            }

            return Ds;
        }

        private double Gamma(double d)
        {
            return d / (2 * d + 1);
        }

        /// <summary>
        /// Equation B_(i,i-1)(u,v) in Eq12, for B_0(u,v) = 1 - GetMPBlending
        /// </summary>
        /// <param name="s"></param>
        /// <param name="d"></param>
        /// <returns></returns>
        private List<double> GetMPBlending(List<double> s, List<double> d)
        {
            List<double> value = new List<double>();
            for (int i = 0; i < d.Count; i++)
            {
                double si = s[i];
                double di = d[i];
                double si_1 = s[IndexWrapper((i - 1), s.Count)];
                double di_1 = d[IndexWrapper((i - 1), d.Count)];

                double vl = (di * HermiteBlend(1 - si_1) * HermiteBlend(di_1) + di_1 * HermiteBlend(si) * HermiteBlend(di)) / (di + di_1);
                value.Add(vl);
            }

            return value;
        }
        private double HermiteBlend(double x)
        {
            return ((Math.Pow((1 - x), 3)) + (3 * (Math.Pow((1 - x), 2)) * x));
        }

        private void SetMesh()
        {
            var trimmedSurf = m_Sph.ToBrep().Split(m_Crvs, DocumentTolerance());
            /*List<Curve> crvs = new List<Curve>();
            for (int i = 0; i < m_Domain.Count; i++)
                crvs.Add(m_Domain[i].ToNurbsCurve());
            var trimmedSurf = m_Sph.ToBrep().Split(crvs, DocumentTolerance());*/

            m_trimmedSphere = new List<Brep>();
            m_trimmedSphere.Add(trimmedSurf[0]);
            m_trimmedSphere.Add(trimmedSurf[1]);

            QuadRemeshParameters meshparams = new QuadRemeshParameters() { TargetQuadCount = 10000 };
            m_DomainMesh = Mesh.QuadRemeshBrep(trimmedSurf[1], meshparams);
        }

        private void SetKatoSurface()
        {
            double u, v;
            List<Point3d> DomainPts = new List<Point3d>();
            List<double> s;
            List<double> d;
            Point3d Pt, P, R;
            double t, tMin, tMax;
            List<double> M;

            for (int i = 0; i < m_Domain.Count; i++)
                DomainPts.Add(m_Domain[i].From);

            for (int i = 0; i < m_DomainMesh.Vertices.Count; i++)
            {
                u = m_DomainMesh.Vertices[i].X;
                v = m_DomainMesh.Vertices[i].Y;
                (s, d) = RadialDistanceFunction(u, v, DomainPts);
                M = GetBlending(d, 1);

                Pt = Point3d.Origin;
                for (int j = 0; j < s.Count; j++)
                {
                    tMin = m_Crvs[j].Domain.Min;
                    tMax = m_Crvs[j].Domain.Max;
                    t = tMin + s[j] * (tMax - tMin);
                    P = m_Crvs[j].PointAt(t);

                    R = P + d[j] * GetRibbonTangent(j, t);
                    Pt = Pt + R * M[j];
                }

                m_DomainMesh.Vertices[i] = (Point3f)Pt;
            }
        }

        private void SetDomain()
        {
            Point3d PtS, PtE;
            m_Domain = new List<Line>();
            for (int i = 0; i < m_Crvs.Count; i++)
            {
                PtS = m_Crvs[i].PointAtStart;
                PtE = m_Crvs[i].PointAtEnd;

                m_Domain.Add(new Line(new Point3d(PtS.X, PtS.Y, 0), new Point3d(PtE.X, PtE.Y, 0)));
            }
        }

        private void SetRibbons()
        {
            m_Sph = new Sphere(new Point3d(0, 0, 0), 2.0);
            m_Sph.Rotate(Math.PI / 2, Vector3d.XAxis);

            List<Point3d> Pts = new List<Point3d>();
            for (int i = 0; i < 6; i++)
                Pts.Add(new Point3d(Math.Cos(i * Math.PI / 3.0), Math.Sin(i * Math.PI / 3.0), 0));

            int j;
            Line Ln;
            Curve[] Curves;
            m_Crvs = new List<Curve>();

            for (int i = 0; i < Pts.Count; i++)
            {
                j = (i + 1) % Pts.Count;
                Ln = new Line(Pts[i], Pts[j]);

                Curves = Rhino.Geometry.Curve.ProjectToBrep(Ln.ToNurbsCurve(), m_Sph.ToBrep(), Vector3d.ZAxis, 0.001);
                if (Curves[0].PointAtEnd.Z > Curves[1].PointAtEnd.Z)
                    m_Crvs.Add(Curves[0]);
                else
                    m_Crvs.Add(Curves[1]);
            }
            SetRibbonTans();
        }

        private void SetRibbonTans()
        {
            m_RibTans = new List<Line>();
            Point3d Pt;
            Vector3d RibTan;
            double tMax, tMin;

            for (int i = 0; i < m_Crvs.Count; i++)
            {
                tMin = m_Crvs[i].Domain.Min;
                tMax = m_Crvs[i].Domain.Max;
                for (double t = tMin; t < tMax + 0.01; t = t + 0.1 * (tMax - tMin))
                {
                    Pt = m_Crvs[i].PointAt(t);
                    RibTan = GetRibbonTangent(i, t);
                    m_RibTans.Add(new Line(Pt, Pt + RibTan));
                }
            }
        }

        private Vector3d GetRibbonTangent(int i, double t)
        {
            Vector3d RibbonTan;

            Vector3d Tan = m_Crvs[i].DerivativeAt(t, 1)[1];
            Point3d Pt = m_Crvs[i].PointAt(t);

            double longitudeRadians, latitudeRadians;
            m_Sph.ClosestParameter(Pt, out longitudeRadians, out latitudeRadians);
            Vector3d Norm = m_Sph.NormalAt(longitudeRadians, latitudeRadians);

            RibbonTan = Vector3d.CrossProduct(Norm, Tan);

            return m_RibScl * RibbonTan;
        }


        /// <summary>
        /// Radial distance functions were suggested by Charrot and Gregory, and these also work for non-regular domains.
        /// </summary>
        /// <param name="u"> u value</param>
        /// <param name="v"> v value</param>
        /// <param name="domainpoint"> domain lines</param>
        /// <param name="si_di"> local parametrization distance from each curve, first item is s value, second item is d value</param>
        private (List<double> si, List<double> di) RadialDistanceFunction(double u, double v, List<Point3d> domainpoint)
        {
            List<double> si = new List<double>();
            List<double> di = new List<double>();
            int i_before, i_after, i_afterafter;

            int n = domainpoint.Count;
            List<Point3d> domainpoints = new List<Point3d>();
            //foreach (var item in domainlines)
            //{
            //    domainpoints.Add(item.From);
            //}
            domainpoints = domainpoint;
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


        /// <summary>
        /// Blending functions
        /// Method >>> 0 - Side Blending, 1 - Special Side Blending, 2 - Corner Blending
        /// </summary>
        /// <param name="d"></param>
        /// <param name="method">0 - Side Blending, 1 - Special Side Blending, 2 - Corner Blending</param>
        /// <returns></returns>
        public List<double> GetBlending(List<double> d, double method)
        {
            List<double> Mus = new List<double>();
            if (method == 0)
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
            if (method == 1)
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
            if (method == 2)
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
                    Numerator = ProductFunction(d, Numerator1);
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
            double Result;
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


        private static int IndexWrapper(int index, int list_count)
        {
            return ((index % list_count) + list_count) % list_count;
        }
        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("2d4e5b59-aa17-49bc-af71-07411b7f57c5"); }
        }
    }
}