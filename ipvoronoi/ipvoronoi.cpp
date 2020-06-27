// includes
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
typedef CGAL::Segment_Delaunay_graph_traits_2<K> Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt> SDG2;
typedef CGAL::Segment_Delaunay_graph_adaptation_traits_2<SDG2> AT;
typedef CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<SDG2> AP;
typedef CGAL::Voronoi_diagram_2<SDG2, AT, AP> VD;
typedef AT::Site_2                    Site_2;
typedef AT::Point_2                   Point_2;
typedef VD::Locate_result             Locate_result;
typedef VD::Vertex_handle             Vertex_handle;
typedef VD::Face_handle               Face_handle;
typedef VD::Halfedge_handle           Halfedge_handle;
typedef VD::Ccb_halfedge_circulator   Ccb_halfedge_circulator;
typedef VD::Bounded_halfedges_iterator BHE_Iter;
typedef VD::Halfedge Halfedge;
typedef VD::Vertex Vertex;
typedef CGAL::Polygon_2<K> Polygon_2;

double distance_to_segment(Point_2 a, Point_2 b, Point_2 x);

int main(int argc, char** argv)
{
  CGAL::set_pretty_mode(std::cout);

  VD vd;
  Site_2 site;
  // Read input file to obtain shape boundary
  FILE *fp = fopen(argv[1], "r");
  std::vector<Point_2> points; points.clear();
  float a, b;

  while (fscanf(fp, "%f %f\n", &a, &b)!=-1) {
      points.push_back(Point_2(a, b));
  }
  fclose(fp);

  // define sites
  for (std::size_t i = 0; i<points.size()-1; i++) {
    site = Site_2::construct_site_2(points[i], points[i+1]);
    vd.insert(site);
  }
  assert( vd.is_valid() );

  // Expensive match between vertex handles and integers. 
  // Then computes edge list as well as edge attributes (up, down) refer to medial axis govenors.
  std::vector<Vertex_handle> vertex_handle_to_int; vertex_handle_to_int.clear();
  std::vector<std::pair<int, int>> edges; edges.clear();
  std::vector<VD::Delaunay_graph::Vertex_handle> ups; ups.clear();
  std::vector<VD::Delaunay_graph::Vertex_handle> downs; downs.clear();
  int uniqueVertCounter = 0;
  BHE_Iter edge_iter = vd.bounded_halfedges_begin();
  int id_count = 0;
  bool found;
  for (;edge_iter != vd.bounded_halfedges_end(); edge_iter++) {
    Halfedge halfedge = *edge_iter;
    Vertex_handle v1p = halfedge.source();
    Vertex_handle v2p = halfedge.target();
    VD::Delaunay_graph::Vertex_handle up = (halfedge.up());
    VD::Delaunay_graph::Vertex_handle down = (halfedge.down());

    int id1 = -1;
    int id2 = -1;
    for (int i = 0; i < vertex_handle_to_int.size(); i++) {
        if (vertex_handle_to_int[i] == v1p) {
            id1 = i;
        }
    }
    if (id1 == -1) {
        id1 = vertex_handle_to_int.size();
        vertex_handle_to_int.push_back(v1p);
    }
    for (int i = 0; i < vertex_handle_to_int.size(); i++) {
        if (vertex_handle_to_int[i] == v2p) {
            id2 = i;
        }
    }
    if (id2 == -1) {
        id2 = vertex_handle_to_int.size();
        vertex_handle_to_int.push_back(v2p);
    }
    edges.push_back(std::make_pair(id1, id2));
    ups.push_back(up);
    downs.push_back(down);
  }

  // Write output details of segment voronoi diagram
  std::ofstream myfile;
  myfile.open("ipvoronoiout.txt");
  for (int i = 0; i < vertex_handle_to_int.size(); i++) {
      myfile << (*(vertex_handle_to_int[i])).point().x() << " " << (*(vertex_handle_to_int[i])).point().y() << std::endl;
  }
  myfile << "edges" << std::endl;
  for (int i = 0; i < edges.size(); i++) {
      myfile << edges[i].first << " " << edges[i].second ;
      myfile << " " << (*(ups[i])).is_point() << " " << (*(downs[i])).is_point() << std::endl;
  }
  myfile << "edgeattributes" << std::endl;
  for (int i = 0; i < ups.size(); i++) {
      auto s1 = (*(ups[i])).site();
      auto s2 = (*(downs[i])).site();
      if (s1.is_point()) {
          myfile << " " << s1.point().x() << " " << s1.point().y() << " ";
      }
      else
      {
          float a1 = s1.segment().source().x();
          float a2 = s1.segment().source().y();
          float a3 = s1.segment().target().x();
          float a4 = s1.segment().target().y();
          myfile << " " << a1 << " " << a2 << " " << a3 << " " << a4 << " " ;
      }

      if (s2.is_point()) {
          myfile << " " << s2.point().x() << " " << s2.point().y() << " ";
      }
      else
      {
          float a1 = s2.segment().source().x();
          float a2 = s2.segment().source().y();
          float a3 = s2.segment().target().x();
          float a4 = s2.segment().target().y();
          myfile << " " << a1 << " " << a2 << " " << a3 << " " << a4 << " ";
      }
      myfile << std::endl;
  }

}
