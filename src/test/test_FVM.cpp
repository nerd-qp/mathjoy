#include "FVM_2dim.hpp"

// // #define DARW_STUFF
#define COMPUTE_STUFF

#include <SFML/Graphics.hpp>

double RandomFloat(double a, double b) {
    const double random = ((double) rand()) / (double) RAND_MAX;
    const double diff = b - a;
    const double r = random * diff;
    return a + r;
}

void benchMark(const std::vector<Vector2<double>>& points,
               const Delaunay<double>& triangulation,
               const std::vector<bool>& on_boundary,
               const utils::function& f,
               const utils::function& target_f) {
  FVM::FVM_solver_2dim test_again(points, triangulation, on_boundary, f);
  test_again.assemble_all();

  // std::cout << test_again.get_stiffness_matrix() << "\n\n";
  // std::cout << test_again.get_right_f_phi() << std::endl;

  auto tmp = test_again.get_solution();
  std::cout << tmp(tmp.size() / 2) << "\n";

  double max_diff = 0;
  size_t max_index = -1;
  double sum = 0;

  auto vertices = test_again.getPoints();

  for ( size_t i = 0; i < vertices.size(); ++ i ) {

    auto results = target_f(vertices[i].x,
                            vertices[i].y);

    std::cout << vertices[i].x << ", "
              << vertices[i].y << " -> "
              << tmp(i) << " vs "
              << results << '\n';

    auto diff = fabs(results - tmp(i));

    sum += max_diff;

    if ( diff > max_diff ) {
      max_index = i;
      max_diff = diff;
    }

    // if ( min_diff > diff ) {
    //   min_diff = diff;
    // }

  }

  std::cout << "Maximum difference is: "
            << max_diff// << ", "
    //<< vertices[max_index].getX() << ", "
    //<< vertices[max_index].getY() << ", "
    //      << max_index
            << '\n';
  std::cout << "Average difference is: "
            << sum / (double)vertices.size() << '\n';
  // std::cout << "Minimum difference is: "
  //           << min_diff << '\n';
}

void good_point_gen(std::vector<Vector2<double>>& points,
                    std::vector<bool>& on_boundary,
                    int numberPoints) {
  const auto n = numberPoints;
  const double deltaX = 1.0 / (n-1);
  const double deltaY = sqrt(3) / (n-1);

  const double x0_1 = 1.0 / (n-1);
  const double x0_2 = 1.5 / (n-1);
  const double y0_1 = 0;
  const double y0_2 = deltaY / 2;

  int stepY = 1.0 / deltaY;

  double xc = x0_1;
  double xc2 = x0_2;
  for ( auto i = 0; i < n - 2; ++ i ) {
    double yc = y0_1;
    double yc2 = y0_2;
    for ( auto j = 0; j < stepY; ++ j ) {
      std::cout << 2 * xc - 1 << ", " << 2 * yc - 1 << '\n';
      std::cout << 2 * xc2 - 1 << ", " << 2 * yc2 - 1 << '\n';
      // std::cout << xc << ", " << yc << '\n';
      // std::cout << xc2 << ", " << yc2 << '\n';
#ifdef COMPUTE_STUFF
      points.push_back(Vector2<double>(2 * xc - 1, 2 * yc - 1));
      points.push_back(Vector2<double>(2 * xc2 - 1, 2 * yc2 - 1));
#endif
#ifdef DARW_STUFF
      points.push_back(Vector2<double>(xc * 800.0, yc * 800.0));
      points.push_back(Vector2<double>(xc2 * 800.0, yc2 * 800.0));
#endif

      yc += deltaY;
      yc2 += deltaY;
    }

    xc += deltaX;
    xc2 += deltaX;
  }

  // add edges points
  xc = 0.0;
  xc2 = 1.0;
  double yc = y0_2;
  for ( auto j = 0u; j < stepY; ++ j ) {
#ifdef COMPUTE_STUFF
    points.push_back(Vector2<double>(2 * xc - 1, 2 * yc - 1));
    points.push_back(Vector2<double>(2 * xc2 - 1, 2 * yc - 1));
#endif
#ifdef DARW_STUFF
    points.push_back(Vector2<double>(xc * 800.0, yc * 800.0));
    points.push_back(Vector2<double>(xc2 * 800.0, yc * 800.0));
#endif

    // points.push_back(Vector2<double>(xc * 800, yc * 800));
    // points.push_back(Vector2<double>(xc2 * 800, yc * 800));
    // points.push_back(Vector2<double>(2 * xc - 1, 2 * yc - 1));
    // points.push_back(Vector2<double>(2 * xc2 - 1, 2 * yc - 1));
    yc += deltaY;
  }

  xc = x0_1;
  xc2 = x0_2;
  yc = 1.0;
  for ( auto i = 0u; i < n; ++ i ) {

#ifdef COMPUTE_STUFF
    points.push_back(Vector2<double>(2 * xc - 1, 2 * yc - 1));
    points.push_back(Vector2<double>(2 * xc2 - 1, 2 * yc - 1));
#endif
#ifdef DARW_STUFF
    points.push_back(Vector2<double>(xc * 800.0, yc * 800.0));
    points.push_back(Vector2<double>(xc2 * 800.0, yc * 800.0));
#endif
    // points.push_back(Vector2<double>(xc * 800, yc * 800));
    // points.push_back(Vector2<double>(xc2 * 800, yc * 800));
    // points.push_back(Vector2<double>(2 * xc - 1, 2 * yc - 1));
    // points.push_back(Vector2<double>(2 * xc2 - 1, 2 * yc - 1));
    xc += deltaX;
    xc2 += deltaX;
  }

  on_boundary.resize(points.size());
  for ( std::size_t i = 0; i < points.size(); ++ i ) {
    const auto& point = points[i];
    // if ( almost_equal(point.x, 0.0) ||
    //      almost_equal(point.x, 1.0) ||
    //      almost_equal(point.y, 0.0) ||
    //      almost_equal(point.y, 1.0) )
    if ( almost_equal(point.x, -1.0) ||
         almost_equal(point.x, 1.0) ||
         almost_equal(point.y, -1.0) ||
         almost_equal(point.y, 1.0) ) {
      //      std::cout << "good\n";
      on_boundary[i] = true;
    }
    else {
      // std::cout << '\n';
      on_boundary[i] = false;
    }
  }

  std::cout << points.size() << ' ' << on_boundary.size() << '\n';
}

void draw(const std::vector<Vector2<double>>& points,
          const Delaunay<double>& triangulation,
          const std::vector<bool>& on_boundary) {
  utils::function f([](double, double) { return 0.0; });
  FVM::FVM_solver_2dim test_run(points, triangulation,
                                  on_boundary,
                                  f);

  const std::vector<Edge<double> >& edges = triangulation.getEdges();

  auto circumCenter_edges = test_run.assemble_all();

  // SFML window
	sf::RenderWindow window(sf::VideoMode(800, 800), "Delaunay triangulation");

	// Transform each points of each vector as a rectangle
	std::vector<sf::RectangleShape*> squares;

	for(const auto p : points) {
		sf::RectangleShape *c1 = new sf::RectangleShape(sf::Vector2f(4, 4));
		c1->setPosition(p.x, p.y);
		squares.push_back(c1);
	}

	// Make the lines
	std::vector<std::array<sf::Vertex, 2> > lines;
	for(const auto &e : edges) {
		lines.push_back({{
			sf::Vertex(sf::Vector2f(e.p1.x + 2, e.p1.y + 2)),
			sf::Vertex(sf::Vector2f(e.p2.x + 2, e.p2.y + 2))
		}});
	}

  std::vector<std::array<sf::Vertex, 2> > circum_lines;
  for(const auto &e : circumCenter_edges) {
		lines.push_back({{
          sf::Vertex(sf::Vector2f(e.p1.x + 2, e.p1.y + 2), sf::Color::Green),
          sf::Vertex(sf::Vector2f(e.p2.x + 2, e.p2.y + 2), sf::Color::Green)
    }});
  }

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		window.clear();

		// Draw the squares
		for(const auto &s : squares) {
			window.draw(*s);
		}

		// Draw the lines
		for(const auto &l : lines) {
			window.draw(l.data(), 2, sf::Lines);
		}

		window.display();
	}

}

int main(int argc, char * argv[])
{
	int numberPoints = 40;
	if (argc==1)
	{
		numberPoints = (int) roundf(RandomFloat(4, numberPoints));
	}
	else if (argc>1)
	{
		numberPoints = atoi(argv[1]);
	}
	std::cout << "Generating " << numberPoints << " good points" << std::endl;

	std::vector<Vector2<double> > points;
  std::vector<bool> on_boundary;
  good_point_gen(points, on_boundary, numberPoints);

	// srand (time(NULL));

	// std::cout << "Generating " << numberPoints << " random points" << std::endl;

	// std::vector<Vector2<double> > points;
	// for(int i = 0; i < numberPoints; ++i) {
	// 	//points.push_back(Vector2<double>(RandomFloat(0, 800), RandomFloat(0, 600)));
  //   points.push_back(Vector2<double>(RandomFloat(0, 1), RandomFloat(0, 1)));
	// }

	Delaunay<double> triangulation;
	const std::vector<Triangle<double> > triangles = triangulation.triangulate(points);
	std::cout << triangles.size() << " triangles generated\n";
	const std::vector<Edge<double> > edges = triangulation.getEdges();

	std::cout << " ========= \n";

  utils::function f([=] (double x, double y) -> double {
                      //return (x*x - 1) * (y*y-1);
                      //return x*x - 1;
                      //return 1;
                      //return M_PI * M_PI * sin(M_PI * x);
                      //return -2;
                      // return - (2 * y * y + 2 * x * x - 4);
                      return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
                      // return sin(M_PI * x) * sin(M_PI * y);
                    });

  utils::function target_f([=] (double x, double y) -> double {
                             //return 1;
                             //return sin(M_PI * x);
                             //return x * x - 1;
                             //return 0;
                             //return (x*x - 1) * (y*y - 1);
                             return sin(M_PI * x) * sin(M_PI * y);
                           });


#ifdef COMPUTE_STUFF
  benchMark(points, triangulation,
            on_boundary, f, target_f);

#endif
#ifdef DARW_STUFF
  draw(points, triangulation, on_boundary);
#endif



}
