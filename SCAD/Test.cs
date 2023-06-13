using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;
using static SCAD.Domain;
using static SCAD.HarmonicClass;

namespace SCAD
{
    public class Test : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Test class.
        /// </summary>
        public Test()
          : base("Transfinite Surface", "Transfinite",
              "Description",
             "SCAD", "Transfinite")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Curves", "C", "The bounding curves", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Domain Type", "D", "Domain Type", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Parametrization", "P", "Parametrization method", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Blending", "B", "Blending Method", GH_ParamAccess.item);
            pManager.AddMeshParameter("Mesh", "M", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Model", "M", "Mesh Model", GH_ParamAccess.item);
            pManager.AddCurveParameter("Domain C", "D", "Domain curves", GH_ParamAccess.list);
            pManager.AddPointParameter("Domain V", "Dv", "Domain vertices", GH_ParamAccess.list);
            pManager.AddMeshParameter("Domain Mesh", "Dm", "Domain Mesh Model", GH_ParamAccess.item);
            pManager.AddPlaneParameter("Ribbon Vector", "Rv", "Ribbon vectors", GH_ParamAccess.list);
            pManager.AddPointParameter("Ribbon center", "Rc", "Ribbon center", GH_ParamAccess.list);
            pManager.AddTextParameter("Ribbon Vector", "Rv", "Ribbon vectors", GH_ParamAccess.tree);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Curve> curves_ = new List<Curve>();
            var ptsDomain = new List<Point3d>();
            var rhinoMesh = new Mesh();
            int d = 0, p = 0, b = 0;
            DA.GetDataList(0, curves_);
            DA.GetData(1, ref d);
            DA.GetData(2, ref p);
            DA.GetData(3, ref b);
            DA.GetData(4, ref rhinoMesh);

            var curves2 = curves_.Select(x => (Curve)x.ToNurbsCurve()).ToList();


            var dm_e = (Domain_Method)d;
            var pm = (Parametrization_Method)p;
            var bm = (Blending_Method)b;
            int resolution = 5;
            #region
            //Domain dm = new Domain();

            //if (dm_e == Domain_Method.Domain_Regular)
            //    dm = new DomainRegular();
            //else if (dm_e == Domain_Method.Domain_Concave)
            //    dm = new DomainConcave();

            //dm.SetSides(curves_);
            //dm.update();
            #endregion

            var surf = new Surface<RibbonCompatible>(dm_e, pm, bm);

            surf.setCurves(curves2);
            surf.setupLoop();
            surf.update();

            #region Visulize domain
            //var dm = surf.GetDomain;
            //var mesh = dm.MeshTopology(resolution);
            //var uvs = dm.Parameters(resolution);
            //var domainMesh = mesh.Getmesh.DuplicateMesh();
            var msh = new Mesh();
            if (dm_e == Domain_Method.Domain_Regular)
            {
                msh = surf.eval(resolution);
                var dm = surf.GetDomain;
                

            }
            else
                msh = surf.eval(rhinoMesh);
            //for (int i = 0; i < uvs.Count; i++)
            //    ptsDomain.Add(new Point3d(uvs[i].X, uvs[i].Y, 0));
            //domainMesh.Vertices.AddVertices(ptsDomain);
            #endregion

            #region
            //var dm = surf.GetDomain;
            //var mesh = dm.MeshTopology(resolution);
            //var uvs = dm.Parameters(resolution);
            //var msh = mesh.Getmesh;
            //var pts = new List<Point3d>();
            //var ptsDomain = new List<Point3d>();
            //var sP = new SurfacePatch(dm, surf.GetRibbons, pm, bm);
            //var domainMesh = msh.DuplicateMesh();
            //for (int i = 0; i < uvs.Count; i++)
            //{
            //    ptsDomain.Add(new Point3d(uvs[i].X, uvs[i].Y, 0));
            //    Point3d katoout = sP.Kato_Suv(uvs[i].X, uvs[i].Y);
            //    pts.Add(katoout);
            //}

            //List<Line> lines = new List<Line>();
            //for (int i = 0; i < sP.vectors.Count; i++)
            //{
            //    lines.Add(new Line(sP.centers[i], sP.vectors[i]));
            //}

            //var sd = dm.Bounds;
            //msh.Vertices.AddVertices(pts);
            //msh.VertexColors.SetColors(Enumerable.Repeat(System.Drawing.Color.Silver, pts.Count).ToArray());
            #endregion

            #region Harmonic map datatree contruction for visualization
            ///
            //var par = surf.GetParametrization;
            //var pointss = domainMesh.Vertices.ToList();
            //DataTree<string> stringtree = new DataTree<string>();
            //for (int i = 0; i < par.HarmonicMapList_si.Count; i++)
            //{
            //    stringtree.Add(writeroutput(par.HarmonicMapList_si[i],pointss), new GH_Path(new int[] { 0, i }));
            //}
            //for (int i = 0; i < par.HarmonicMapList_di.Count; i++)
            //{
            //    stringtree.Add(writeroutput(par.HarmonicMapList_di[i],pointss), new GH_Path(new int[] { 1, i }));
            //}
            ///
            #endregion

            DA.SetData(0, msh);
            //DA.SetDataList(2, dm.Vertices);
            //DA.SetData(3, domainMesh);
            //DA.SetDataList(4, sP.planes);
            //DA.SetDataList(5, pointss);
            //DA.SetDataTree(6, stringtree);

        }

        //private string writeroutput(HarmonicMap map, List<Point3f> pointss)
        //{
        //    string f = "";
        //    for (int i = 0; i < pointss.Count; i++)
        //    {
        //        Point3d point2 = pointss[i];
        //        var success = harmonic_eval(map, point2, out double result);
        //        if (success)
        //        {
        //            //f.Write("{0:N6},{0:N6},{0:N6}\n", u, v, result);
        //            f += (point2.X.ToString("N6") + "," + point2.Y.ToString("N6") + "," + result.ToString("N6") + '\n');
        //            //f.Write(u.ToString() + ',' + v.ToString() + ',' + result.ToString() + '\n');
        //        }
        //    }
        //    return f;
        //}

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
            get { return new Guid("5587F090-35B5-45DE-8AE8-64D81549853E"); }
        }
    }
}