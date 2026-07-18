#include <fstream>
#include <vector>
#include <map>

struct Point {
    double x, y;
};

using Polygon = std::vector<Point>;

std::map<int, Polygon> read_regions(const std::string& filename)
{
    std::ifstream in(filename);
    std::map<int, Polygon> regions;

    int id;
    double x, y, a, b;

    while (in >> id >> x >> y >> a >> b) {
        regions[id].push_back({x, y});
    }

    return regions;
}

bool point_in_polygon(double x, double y, const Polygon& poly)
{
    bool inside = false;
    size_t n = poly.size();

    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        double xi = poly[i].x, yi = poly[i].y;
        double xj = poly[j].x, yj = poly[j].y;

        bool intersect =
            ((yi > y) != (yj > y)) &&
            (x < (xj - xi) * (y - yi) / (yj - yi) + xi);

        if (intersect)
            inside = !inside;
    }
    return inside;
}

void loop_over_region(const Polygon& poly, int id, double dx)
{
    double xmin = poly[0].x, xmax = poly[0].x;
    double ymin = poly[0].y, ymax = poly[0].y;
    std::ofstream out("points.dat");

    for (const auto& p : poly) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
    }

    for (double x = xmin; x <= xmax; x += dx) {
        for (double y = ymin; y <= ymax; y += dx) {
            if (point_in_polygon(x, y, poly)) {
                out << x << " " << y << " " << id << "\n";
            }
        }
    }
    out.close();
}

int main()
{
    auto regions = read_regions("layout_7f_3.outline.lbad");

    double dx = 0.01;   // spacing between points

    for (const auto& [id, poly] : regions) {
        loop_over_region(poly, id, dx);
    }

    return 0;
}
