#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/hierarchy_simplify_point_set.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Timer.h>

#include <fstream>
#include <list>
#include <cassert>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Alpha_shape_vertex_base_3<K>               Vb;
typedef CGAL::Alpha_shape_cell_base_3<K>                 Fb;
typedef CGAL::Triangulation_data_structure_3<Vb, Fb>      Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>  Delaunay;
typedef CGAL::Alpha_shape_3<Delaunay>                    Alpha_shape_3;
typedef K::Point_3                                       Point;
typedef Alpha_shape_3::Alpha_iterator                    Alpha_iterator;
typedef Alpha_shape_3::NT                                NT;

void writeOFF(Alpha_shape_3 &as, std::string asfile)
{
	/// collect all regular facets
	std::vector<Alpha_shape_3::Facet> facets;
	as.get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);
	std::cout << facets.size() << " triangles" << std::endl;

	std::stringstream pts;
	std::stringstream ind;

	std::size_t nbf = facets.size();
	for (std::size_t i = 0; i < nbf; ++i)
	{
		//To have a consistent orientation of the facet, always consider an exterior cell
		if (as.classify(facets[i].first) != Alpha_shape_3::EXTERIOR)
			facets[i] = as.mirror_facet(facets[i]);
		CGAL_assertion(as.classify(facets[i].first) == Alpha_shape_3::EXTERIOR);

		int indices[3] = {
			(facets[i].second + 1) % 4,
			(facets[i].second + 2) % 4,
			(facets[i].second + 3) % 4,
		};

		/// according to the encoding of vertex indices, this is needed to get
		/// a consistent orienation
		if (facets[i].second % 2 == 0) std::swap(indices[0], indices[1]);

		pts <<
			facets[i].first->vertex(indices[0])->point() << "\n" <<
			facets[i].first->vertex(indices[1])->point() << "\n" <<
			facets[i].first->vertex(indices[2])->point() << "\n";
		ind << "3 " << 3 * i << " " << 3 * i + 1 << " " << 3 * i + 2 << "\n";
	}

	std::ofstream fout(asfile);
	fout << "OFF\n" << 3 * nbf << " " << nbf << " 0\n";
	fout << pts.str();
	fout << ind.str();
	fout.close();
}

int main()
{
	for (int i = 1; i < 3; ++i)
	{
		for (int mat = 1; mat < 2; ++mat)
		{
			std::string data_path = "F:/project/def/def_" + std::to_string(i+1) + "_" + std::to_string(mat+1) + ".xyz";

			std::vector<Point> points;
			std::ifstream stream(data_path);
			if (!stream ||
				!CGAL::read_xyz_points(stream, std::back_inserter(points)))
			{
				continue;
			}
			std::cout << "Read " << points.size() << " point(s)" << std::endl;
			CGAL::Timer task_timer; task_timer.start();

			unsigned int n = points.size();
			unsigned int para_n = 10000;
			if (n < 1000000) para_n = 1000;
			if (n < 100000) para_n = 100;
			if (n < 10000) para_n = 10;

			/// simplification by clustering using erase-remove idiom
			points.erase(CGAL::hierarchy_simplify_point_set(points,
						CGAL::parameters::size(para_n). // Max cluster size
						maximum_variation(0.001)), // Max surface variation
						points.end());

			std::size_t memory = CGAL::Memory_sizer().virtual_size();
			std::cout << points.size() << " point(s) kept, computed in "
				<< task_timer.time() << " seconds, "
				<< (memory >> 20) << " Mib allocated." << std::endl;

			Delaunay dt;
			Point p;
			
			for (; n > 0; n--) {
				dt.insert(points[n]);
			}
			std::cout << "Delaunay computed." << std::endl;

			/// compute alpha shape
			Alpha_shape_3 as(dt);
			std::cout << "Alpha shape computed in REGULARIZED mode by default." << std::endl;

			/// 1. find smallest alpha values to get a solid, not always connected
			// Alpha_shape_3::NT alpha_solid = as.find_alpha_solid();
			// std::cout << "Smallest alpha value to get a solid through data points is " << alpha_solid << std::endl;
			// as.set_alpha(alpha_solid);

			/// 2. optimal alpha values to get connected componet
			Alpha_iterator opt = as.find_optimal_alpha(1);
			std::cout << "Optimal alpha value to get one connected component is " << *opt << std::endl;
			//Alpha_shape_3::NT alpha_final = (*opt + alpha_solid) / 2;
			as.set_alpha(*opt);

			std::string asfile = "F:/project/mesh/" + std::to_string(i + 1) + "_" + std::to_string(mat + 1) + ".off";
			writeOFF(as, asfile);
		}
	}

	system("pause");
	return 0;
}