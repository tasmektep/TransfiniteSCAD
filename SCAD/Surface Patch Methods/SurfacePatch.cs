using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper;
using Rhino.Geometry;

namespace SCAD
{
    class SurfacePatch
    {

        ////
        ///Blending function bana di değeri girdisi verildiğinde n (yani curve sayısı) kadar value değeri vermeli.

        /// <summary>
        /// kato
        /// curve üzerinde ribbonlara ihtiyaç var, eğer curve içinde gömülü değilse
        /// </summary>
        /// <param name="u"> u value </param>
        /// <param name="v"> v value </param>
        /// <param name="curves"> boundary curves </param>
        /// <returns></returns>
        //private Point3d Kato_Suv(double u, double v, List<Curve> curves)
        //{
        //    new Domain(curves,out List<Curve> domaincurve);
        //    //new RadialDistanceFunction(u, v, domaincurve, out List<(double, double)> si_di);
        //    var prm = new Parametrization(Parametrization_Method.Harmonic);


        //    List<(double, double)> si_di = prm.GetPoint(u, v);
        //    ///Alper sence bu parametrization seçimi nerede yapılmalı parametrization tarafında yapılıp oradan sadece tek bir si_di listi mi verilmeli.

        //    //List<(double, double)> si_di = ComputeDistance2(u, v);
        //    //List<(double, double)> si_di = ComputeDistance3(u, v);
        //    //List<(double, double)> si_di = MVC(u, v);
        //    //List<(double, double)> si_di = HarmonicMapCall(u, v, m_Curves.Count);

        //    //List<(double, double)> si_di = HarmonicMapCall(u, v, m_CurvesSerhat.Count);//curve part
        //    //m_uvSerhat.Add(new Point3d(u, v, 0));
        //    //int j;
        //    //Line Ln;
        //    //Point3d Pt;
        //    //double s, d;
        //    /* for (int i = 0; i < m_DomainPolygon.Count; i++)
        //     {
        //         j = (i+1)% m_DomainPolygon.Count;
        //         Ln = new Line(m_DomainPolygon[i], m_DomainPolygon[j]);
        //         Pt = Ln.ClosestPoint(new Point3d(u, v, 0), true);
        //         s = Pt.DistanceTo(Ln.From) / Ln.Length;
        //         d = Pt.DistanceTo(new Point3d(u, v, 0));

        //         si_di.Add((s, d));
        //     }*/

        //    double s;
        //    List<double> d_i = new List<double>();
        //    for (int i = 0; i < si_di.Count; i++)
        //    {
        //        d_i.Add(si_di[i].Item2);
        //    }

        //    Vector3d T = new Vector3d();
        //    Point3d r_sum = new Point3d();
        //    //List<double> Value = ComputeSPsideBLFunction1(d_i);
        //    //BlendingFunctions value = new BlendingFunctions(0, BlendingFunctions.Blending_Methods.Special_Side_Blending); ** Arashdan gelen data*** list value olarak gelse güzel olur.
        //    int j;
        //    for (int i = 0; i < curves.Count; i++)
        //    {
        //        //double domainlengt = (m_DomainPolygon[IndexWrapper(i, m_Curves.Count)] - m_DomainPolygon[IndexWrapper(i + 1, m_Curves.Count)]).Length;
        //        double curveparameter = (curves[i].Domain.Length * si_di[i].Item1) / curves[i].GetLength();

        //        //Vector3d crossproduct = Vector3d.CrossProduct(m_Curves[i].TangentAt(curveparameter), m_Curves[i].CurvatureAt(curveparameter));
        //        //Point3d r = m_Curves[i].PointAt(curveparameter) + (si_di[i].Item2 * crossproduct);
        //        //s = m_Curves[i].Domain.Min + si_di[i].Item1 * (m_Curves[i].Domain.Max - m_Curves[i].Domain.Min);
        //        //Vector3d crossproduct = Vector3d.CrossProduct(m_Curves[i].TangentAt(s), m_Curves[i].CurvatureAt(s));
        //        //Vector3d crossproduct = m_Curves[i].CurvatureAt(s);
        //        //Point3d r = m_Curves[i].PointAt(s) + (si_di[i].Item2 * crossproduct);

        //        //Point3d r = m_Curves[i].PointAt(curveparameter) + (si_di[i].Item2 * m_Curves[i].TangentAt(curveparameter));
        //        //   s = m_Curves[i].Domain.Min + si_di[i].Item1 * (m_Curves[i].Domain.Max - m_Curves[i].Domain.Min);
        //        //   Point3d r = m_Curves[i].PointAt(s) + (si_di[i].Item2 * m_Curves[i].TangentAt(s));
        //        //double blendingfunctionvalue = ComputeSPsideBLFunction(d_i, i);

        //        s = curves[i].Domain.Min + si_di[i].Item1 * (curves[i].Domain.Max - curves[i].Domain.Min);
        //        j = (i + 1) % curves.Count;
        //        //T = m_Ts[i] + si_di[i].Item1 * (m_Ts[j] - m_Ts[i]);
        //        //if (m_TssSerhat[i].Count > 2)
        //        //{
        //        //    ;
        //        //}
        //        //else T = m_TssSerhat[i][0] + si_di[i].Item1 * (m_TssSerhat[i][1] - m_TssSerhat[i][0]);

        //        T = m_TssSerhat[i][0] + (si_di[i].Item1 / (curves[i].Domain.Max - curves[i].Domain.Min)) * (m_TssSerhat[i][1] - m_TssSerhat[i][0]);


        //        Point3d r = curves[i].PointAt(s) + (si_di[i].Item2 * T);

        //        r_sum += r * Value[i];
        //    }

        //    return r_sum;
        //}
    }
}
