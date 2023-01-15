using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SCAD
{
    public static class Extensions
    {
        public static int IndexWrapper(int index, int list_count)
        {
            return ((index % list_count) + list_count) % list_count;
        }
    }
}
