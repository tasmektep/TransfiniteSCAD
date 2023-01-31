using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace SCAD
{
    public class Domain
    {


        private 
            List<Curve> m_BndCurves;
            List<Curve> m_DomainCurves;

        /// <summary>
        /// Gets boundary curves to create Domain Curves
        /// returns: The Domain Curves
        /// </summary>
        /// <param name="BndCurves"> Boundary Curves </param>
        /// <param name="DomainCurves"> Domain Curves </param>
        public Domain(List<Curve> BndCurves, out List<Curve> DomainCurves)
        {
            m_BndCurves = BndCurves;
            ComputeDomainPolygon();
            DomainCurves = m_DomainCurves;
        }

        private void ComputeDomainPolygon()
        {
            double TotL = 0, L;
            List<double> Ls = new List<double>();
            for (int i = 0; i < m_BndCurves.Count; i++)
            {
                L = m_BndCurves[i].GetLength();
                TotL = TotL + L;
                Ls.Add(L);
            }


            List<double> Alphas = new List<double>();
            Alphas.Add(0.0);
            for (int i = 1; i < m_BndCurves.Count; i++)
            {
                L = 0;
                for (int j = 0; j <= i - 1; j++)
                    L = L + Ls[j];
                Alphas.Add(2 * Math.PI * L / TotL);
            }

            m_DomainCurves = new List<Curve>();
            Point3d Pt1, Pt2;
            int k;
            Line Ln;
            for (int i = 0; i < Alphas.Count; i++)
            {
                k = (i + 1) % Alphas.Count;
                Pt1 = new Point3d(Math.Cos(Alphas[i]), Math.Sin(Alphas[i]), 0);
                Pt2 = new Point3d(Math.Cos(Alphas[k]), Math.Sin(Alphas[k]), 0);

                Ln = new Line(Pt1, Pt2);
                m_DomainCurves.Add(Ln.ToNurbsCurve());
            }
        }
    }

}
