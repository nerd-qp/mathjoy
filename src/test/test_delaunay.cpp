#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <array>

#include <SFML/Graphics.hpp>

#include "delaunay-triangulation/vector2.h"
#include "delaunay-triangulation/triangle.h"
#include "delaunay-triangulation/delaunay.h"

double RandomFloat(double a, double b) {
    const double random = ((double) rand()) / (double) RAND_MAX;
    const double diff = b - a;
    const double r = random * diff;
    return a + r;
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
	srand (time(NULL));

	// std::cout << "Generating " << numberPoints << " random points" << std::endl;
  std::cout << "Generating " << numberPoints << " good points" << std::endl;

	std::vector<Vector2<double> > points;
  const auto n = numberPoints;
  const double deltaX = 1.0 / (n-1);
  const double deltaY = sqrt(3) / (n-1);

  const double x0_1 = 0;
  const double x0_2 = 0.5 / (n-1);
  const double y0_1 = 0;
  const double y0_2 = deltaY / 2;

  int stepY = 1.0 / deltaY;

  double xc = x0_1;
  double xc2 = x0_2;
  for ( auto i = 0; i < n; ++ i ) {
    double yc = y0_1;
    double yc2 = y0_2;
    for ( auto j = 0; j < stepY; ++ j ) {
      std::cout << xc << ", " << yc << '\n';
      std::cout << xc2 << ", " << yc2 << '\n';
      points.push_back(Vector2<double>(xc * 800.0, yc * 600.0));
      points.push_back(Vector2<double>(xc2 * 800.0, yc2 * 600.0));
      yc += deltaY;
      yc2 += deltaY;
    }
    xc += deltaX;
    xc2 += deltaX;
  }

	for(int i = 0; i < numberPoints; ++i) {
		//points.push_back(Vector2<double>(RandomFloat(0, 800), RandomFloat(0, 600)));
	}

	Delaunay<double> triangulation;
	const std::vector<Triangle<double> > triangles = triangulation.triangulate(points);
	std::cout << triangles.size() << " triangles generated\n";
	const std::vector<Edge<double> > edges = triangulation.getEdges();

	std::cout << " ========= ";

	// std::cout << "\nPoints : " << points.size() << std::endl;
	// for(const auto &p : points)
	// 	std::cout << p << std::endl;

	// std::cout << "\nTriangles : " << triangles.size() << std::endl;
	// for(const auto &t : triangles)
	// 	std::cout << t << std::endl;

	// std::cout << "\nEdges : " << edges.size() << std::endl;
	// for(const auto &e : edges)
	// 	std::cout << e << std::endl;

	// SFML window
	sf::RenderWindow window(sf::VideoMode(800, 600), "Delaunay triangulation");

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

	return 0;
}
