using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper.Kernel;
using Rhino.Geometry;
using static SCAD.Domain;

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
            pManager.AddVectorParameter("Ribbon Vector", "Rv", "Ribbon vectors", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Curve> curves_ = new List<Curve>();
            int d = 0, p = 0, b = 0;
            DA.GetDataList(0, curves_);
            DA.GetData(1, ref d);
            DA.GetData(2, ref p);
            DA.GetData(3, ref b);

            var dm_e = (Domain_Method)d;
            var pm = (Parametrization_Method)p;
            var bm = (Blending_Method)b;
            int resolution = 10;
            Domain dm = new Domain();

            if (dm_e == Domain_Method.Domain_Regular)
                dm = new DomainRegular();
            else if (dm_e == Domain_Method.Domain_Curved)
                dm = new DomainCurved();

            dm.SetSides(curves_);
            dm.update();

            var surf = new SCAD.Surface<DomainCurved, Parametrization, RibbonCompatible>();
            double scaling = 20.0;
            double ribbon_length = 0.25;

            surf.setCurves(curves_);
            surf.setupLoop();
            surf.update();

            var mesh = dm.MeshTopology(resolution);
            var uvs = dm.Parameters(resolution);
            var msh = mesh.Getmesh;
            var pts = new List<Point3d>();
            var ptsDomain = new List<Point3d>();
            var sP = new SurfacePatch(dm, surf.GetRibbons, pm, bm);
            var domainMesh = msh.DuplicateMesh();
            for (int i = 0; i < uvs.Count; i++)
            {
                ptsDomain.Add(new Point3d(uvs[i].X, uvs[i].Y, 0));
                Point3d katoout = sP.Kato_Suv(uvs[i].X, uvs[i].Y);
                pts.Add(katoout);
            }

            List<Vector3d> vecList = new List<Vector3d>();
            List<Point3d> ptList = new List<Point3d>();
            for (int i = 0; i < surf.n; ++i)
            {
                for (int j = 0; j <= resolution; ++j)
                {
                    double u = (double)j / resolution;
                    ptList.Add(surf.ribbon(i).curve().PointAt(u));
                    vecList.Add((Vector3d)surf.ribbon(i).eval(new Point2d(u, ribbon_length)));
                }
            }





            //List<Line> lines = new List<Line>();
            //for (int i = 0; i < sP.vectors.Count; i++)
            //{
            //    lines.Add(new Line(sP.centers[i], sP.vectors[i]));
            //}

            //var sd = dm.Bounds;
            msh.Vertices.AddVertices(pts);
            domainMesh.Vertices.AddVertices(ptsDomain);
            msh.VertexColors.SetColors(Enumerable.Repeat(System.Drawing.Color.Silver, pts.Count).ToArray());
            DA.SetData(0, msh);
            DA.SetDataList(2, dm.Vertices);
            DA.SetData(3, domainMesh);
            DA.SetDataList(4, sP.planes);
            //DA.SetDataList(5, ptList);
            //DA.SetDataList(6, vecList);

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
            get { return new Guid("5587F090-35B5-45DE-8AE8-64D81549853E"); }
        }
    }
}