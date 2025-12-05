#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Modifier_base.h> 

#include <vector>
#include <list>
#include <iostream>
#include <string>
#include <algorithm>
#include <map>
#include <set>
#include <memory>
#include <cmath>
#include <random>
#include <chrono>
#include <sstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef K::Point_3 Point_3;
typedef K::Point_2 Point_2;
typedef K::Vector_3 Vector_3;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::HalfedgeDS HalfedgeDS;

// =========================================================
// PART 1: GENERATOR
// =========================================================

class NoiseModifier : public CGAL::Modifier_base<HalfedgeDS> {
    double phase_x, phase_y, phase_z;     // Зсув для форми
    double phase_x2, phase_y2, phase_z2;  // Зсув для шипів
public:
    NoiseModifier() {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine rng(seed);
        std::uniform_real_distribution<double> dist(0.0, 100.0);
        
        phase_x = dist(rng); phase_y = dist(rng); phase_z = dist(rng);
        phase_x2 = dist(rng); phase_y2 = dist(rng); phase_z2 = dist(rng);
    }

    void operator()(HalfedgeDS& hds) {
        for (auto v = hds.vertices_begin(); v != hds.vertices_end(); ++v) {
            double x = CGAL::to_double(v->point().x());
            double y = CGAL::to_double(v->point().y());
            double z = CGAL::to_double(v->point().z());
            
            double len = std::sqrt(x*x + y*y + z*z);
            if (len < 1e-9) continue;
            x /= len; y /= len; z /= len;
            
            double shape = sin(2.0*x + phase_x) * cos(2.0*y + phase_y) * sin(2.0*z + phase_z);
            
            double spikes = sin(5.0*x + phase_x2) + cos(5.0*y + phase_y2) + sin(5.0*z + phase_z2);
            
            double displacement = 1.0 + 0.4 * shape + 0.3 * spikes;
            
            displacement = std::max(0.2, displacement);

            v->point() = Point_3(x * displacement, y * displacement, z * displacement);
        }
    }
};

void generate_star_mesh(int num_points, Polyhedron& result_mesh) {
    std::vector<Point_3> points;
    CGAL::Random_points_on_sphere_3<Point_3> generator(1.0);
    CGAL::cpp11::copy_n(generator, num_points, std::back_inserter(points));

    Polyhedron primal;
    CGAL::convex_hull_3(points.begin(), points.end(), primal);

    // Dual Mesh Construction
    std::vector<Point_3> dual_vertices;
    std::vector<std::vector<int>> dual_faces;
    std::map<Polyhedron::Facet_const_handle, int> facet_to_idx;

    int idx = 0;
    for (auto f = primal.facets_begin(); f != primal.facets_end(); ++f) {
        auto h = f->halfedge();
        Point_3 p0 = h->vertex()->point();
        Point_3 p1 = h->next()->vertex()->point();
        Point_3 p2 = h->next()->next()->vertex()->point();
        
        double cx = (CGAL::to_double(p0.x()) + CGAL::to_double(p1.x()) + CGAL::to_double(p2.x())) / 3.0;
        double cy = (CGAL::to_double(p0.y()) + CGAL::to_double(p1.y()) + CGAL::to_double(p2.y())) / 3.0;
        double cz = (CGAL::to_double(p0.z()) + CGAL::to_double(p1.z()) + CGAL::to_double(p2.z())) / 3.0;
        
        dual_vertices.emplace_back(cx, cy, cz);
        facet_to_idx[f] = idx++;
    }

    for (auto v = primal.vertices_begin(); v != primal.vertices_end(); ++v) {
        std::vector<int> face_indices;
        auto h_start = v->vertex_begin();
        auto h_circ = h_start;
        do {
            if (!h_circ->is_border()) face_indices.push_back(facet_to_idx[h_circ->facet()]);
            h_circ++;
        } while (h_circ != h_start);

        if (face_indices.size() >= 3) {
            std::reverse(face_indices.begin(), face_indices.end());
            dual_faces.push_back(face_indices);
        }
    }

    std::stringstream ss;
    ss << "OFF\n" << dual_vertices.size() << " " << dual_faces.size() << " 0\n";
    for(const auto& p : dual_vertices) ss << p << "\n";
    for(const auto& f : dual_faces) {
        ss << f.size();
        for(int i : f) ss << " " << i;
        ss << "\n";
    }
    ss >> result_mesh;

    NoiseModifier modifier;
    result_mesh.delegate(modifier);
}

// =========================================================
// PART 2: SWEEP-LINE
// =========================================================

struct DcelHalfEdge;
struct DcelVertex {
    int id; Point_2 p; DcelHalfEdge* incident_edge = nullptr; Halfedge_handle original_h; 
};
struct DcelHalfEdge {
    int id; DcelVertex* origin = nullptr; DcelHalfEdge* next = nullptr; DcelHalfEdge* prev = nullptr;
    mutable DcelVertex* helper = nullptr; 
};

class DCEL {
public:
    std::list<DcelVertex> vertices;
    std::list<DcelHalfEdge> edges;

    DcelVertex* createVertex(int id, Point_2 p, Halfedge_handle h) {
        vertices.push_back(DcelVertex());
        vertices.back().id = id; 
        vertices.back().p = p; 
        vertices.back().original_h = h;
        return &vertices.back();
    }
    DcelHalfEdge* createEdge(int id) {
        edges.push_back(DcelHalfEdge()); 
        edges.back().id = id; 
        return &edges.back();
    }
};

enum VertexType { START, END, SPLIT, MERGE, REGULAR_LEFT, REGULAR_RIGHT };
struct Event {
    DcelVertex* v;
    bool operator<(const Event& other) const {
        double y1 = CGAL::to_double(v->p.y()); double y2 = CGAL::to_double(other.v->p.y());
        if (std::abs(y1 - y2) > 1e-9) return y1 > y2; 
        return CGAL::to_double(v->p.x()) < CGAL::to_double(other.v->p.x()); 
    }
};

VertexType getType(DcelVertex* v, DcelVertex* prev, DcelVertex* next) {
    double vy = CGAL::to_double(v->p.y()); double py = CGAL::to_double(prev->p.y()); double ny = CGAL::to_double(next->p.y());
    auto is_below = [&](double y_test, double y_ref, double x_test, double x_ref) {
        return (y_test < y_ref) || (std::abs(y_test - y_ref) < 1e-9 && x_test > x_ref);
    };
    bool prev_is_above = !is_below(py, vy, CGAL::to_double(prev->p.x()), CGAL::to_double(v->p.x()));
    bool next_is_above = !is_below(ny, vy, CGAL::to_double(next->p.x()), CGAL::to_double(v->p.x()));
    if (prev_is_above && next_is_above) {
        if (CGAL::orientation(prev->p, v->p, next->p) == CGAL::RIGHT_TURN) return MERGE; return END;
    }
    if (!prev_is_above && !next_is_above) {
        if (CGAL::orientation(prev->p, v->p, next->p) == CGAL::RIGHT_TURN) return SPLIT; return START;
    }
    if (prev_is_above) return REGULAR_RIGHT; return REGULAR_LEFT;
}

struct StatusEdge {
    DcelHalfEdge* edge;
    double getX(double currentY) const {
        Point_2 p1 = edge->origin->p; Point_2 p2 = edge->next->origin->p;
        double y1 = CGAL::to_double(p1.y()), y2 = CGAL::to_double(p2.y());
        double x1 = CGAL::to_double(p1.x()), x2 = CGAL::to_double(p2.x());
        if (std::abs(y1 - y2) < 1e-9) return std::max(x1, x2);
        double t = (currentY - y1) / (y2 - y1);
        return x1 + t * (x2 - x1);
    }
};
double SWEEP_Y = 0.0;
struct StatusComp {
    bool operator()(const StatusEdge* a, const StatusEdge* b) const {
        double xa = a->getX(SWEEP_Y); double xb = b->getX(SWEEP_Y);
        if (std::abs(xa - xb) < 1e-9) return a->edge->id < b->edge->id;
        return xa < xb;
    }
};

void triangulate_facet_fan(Facet_handle f, Polyhedron& mesh) {
    Halfedge_handle h = f->halfedge(); int max_iter = 1000;
    while (h->next()->next()->next() != h && max_iter-- > 0) { mesh.split_facet(h, h->next()->next()); }
}

bool triangulate_facet_sweepline(Facet_handle f, Polyhedron& mesh) {
    if (f->is_triangle()) return true;
    DCEL dcel;
    Halfedge_handle h_start = f->halfedge(), h = h_start;
    double nx = 0, ny = 0, nz = 0;
    do {
        Point_3 p1 = h->vertex()->point(), p2 = h->next()->vertex()->point();
        double x1=CGAL::to_double(p1.x()), y1=CGAL::to_double(p1.y()), z1=CGAL::to_double(p1.z());
        double x2=CGAL::to_double(p2.x()), y2=CGAL::to_double(p2.y()), z2=CGAL::to_double(p2.z());
        nx += (y1-y2)*(z1+z2); ny += (z1-z2)*(x1+x2); nz += (x1-x2)*(y1+y2);
        h = h->next();
    } while (h != h_start);
    int drop = (std::abs(nx)>std::abs(ny) && std::abs(nx)>std::abs(nz)) ? 0 : (std::abs(ny)>std::abs(nz) ? 1 : 2);

    std::vector<DcelVertex*> poly_vertices; int idx = 0; h = h_start;
    do {
        Point_3 p3 = h->vertex()->point(); Point_2 p2;
        if (drop == 0) p2 = Point_2(p3.y(), p3.z()); else if (drop == 1) p2 = Point_2(p3.x(), p3.z()); else p2 = Point_2(p3.x(), p3.y());
        poly_vertices.push_back(dcel.createVertex(idx++, p2, h)); h = h->next();
    } while (h != h_start);

    int N = poly_vertices.size(); if (N <= 3) return true; 
    for (int i = 0; i < N; ++i) {
        DcelHalfEdge* e = dcel.createEdge(i); e->origin = poly_vertices[i]; poly_vertices[i]->incident_edge = e;
    }
    auto it_e = dcel.edges.begin();
    for(int i=0; i<N; ++i) {
        auto prev_it = (i == 0) ? std::prev(dcel.edges.end()) : std::prev(it_e);
        auto next_it = (i == N-1) ? dcel.edges.begin() : std::next(it_e);
        it_e->prev = &(*prev_it); it_e->next = &(*next_it); it_e++;
    }

    std::vector<Event> events; for(auto v : poly_vertices) events.push_back({v});
    std::sort(events.begin(), events.end());

    std::set<StatusEdge*, StatusComp> status;
    std::map<int, StatusEdge*> edge_map;
    std::vector<std::pair<Halfedge_handle, Halfedge_handle>> diagonals;
    std::list<StatusEdge> status_storage; 

    for (const auto& ev : events) {
        DcelVertex* v = ev.v; SWEEP_Y = CGAL::to_double(v->p.y());
        DcelVertex* prev = v->incident_edge->prev->origin; DcelVertex* next = v->incident_edge->next->origin;
        VertexType type = getType(v, prev, next);
        DcelHalfEdge* e_v = v->incident_edge; DcelHalfEdge* e_prev = v->incident_edge->prev; 
        try {
            if (type == START) {
                status_storage.push_back({e_v}); StatusEdge* st = &status_storage.back();
                st->edge->helper = v; status.insert(st); edge_map[e_v->id] = st;
            } else if (type == END) {
                if (edge_map.count(e_prev->id)) {
                    StatusEdge* st = edge_map[e_prev->id];
                    if (st->edge->helper && getType(st->edge->helper, st->edge->helper->incident_edge->prev->origin, st->edge->helper->incident_edge->next->origin) == MERGE)
                        diagonals.push_back({v->original_h, st->edge->helper->original_h});
                    status.erase(st); edge_map.erase(e_prev->id);
                } else return false;
            } else if (type == SPLIT) {
                StatusEdge probe_obj{e_v}; auto it = status.lower_bound(&probe_obj);
                if (it != status.begin()) {
                    --it; StatusEdge* left = *it;
                    diagonals.push_back({v->original_h, left->edge->helper->original_h});
                    left->edge->helper = v;
                }
                status_storage.push_back({e_v}); StatusEdge* st = &status_storage.back();
                st->edge->helper = v; status.insert(st); edge_map[e_v->id] = st;
            } else if (type == MERGE) {
                if (edge_map.count(e_prev->id)) {
                    StatusEdge* st = edge_map[e_prev->id];
                    if (st->edge->helper && getType(st->edge->helper, st->edge->helper->incident_edge->prev->origin, st->edge->helper->incident_edge->next->origin) == MERGE)
                        diagonals.push_back({v->original_h, st->edge->helper->original_h});
                    status.erase(st); edge_map.erase(e_prev->id);
                } else return false;
                StatusEdge probe_obj{e_v}; auto it = status.lower_bound(&probe_obj);
                if (it != status.begin()) {
                    --it; StatusEdge* left = *it;
                    if (left->edge->helper && getType(left->edge->helper, left->edge->helper->incident_edge->prev->origin, left->edge->helper->incident_edge->next->origin) == MERGE)
                         diagonals.push_back({v->original_h, left->edge->helper->original_h});
                    left->edge->helper = v;
                }
            } else if (type == REGULAR_LEFT) {
                if (edge_map.count(e_prev->id)) {
                    StatusEdge* st = edge_map[e_prev->id];
                    if (st->edge->helper && getType(st->edge->helper, st->edge->helper->incident_edge->prev->origin, st->edge->helper->incident_edge->next->origin) == MERGE)
                        diagonals.push_back({v->original_h, st->edge->helper->original_h});
                    status.erase(st); edge_map.erase(e_prev->id);
                }
                status_storage.push_back({e_v}); StatusEdge* st = &status_storage.back();
                st->edge->helper = v; status.insert(st); edge_map[e_v->id] = st;
            } else if (type == REGULAR_RIGHT) {
                StatusEdge probe_obj{e_v}; auto it = status.lower_bound(&probe_obj);
                if (it != status.begin()) {
                    --it; StatusEdge* left = *it;
                    if (left->edge->helper && getType(left->edge->helper, left->edge->helper->incident_edge->prev->origin, left->edge->helper->incident_edge->next->origin) == MERGE)
                        diagonals.push_back({v->original_h, left->edge->helper->original_h});
                    left->edge->helper = v;
                }
            }
        } catch (...) { return false; } 
    }
    try {
        for (auto& pair : diagonals) 
            if (pair.first->facet() == pair.second->facet()) mesh.split_facet(pair.first, pair.second);
    } catch(...) { return false; }
    return true;
}

void triangulate_monotone_cleanup(Facet_handle f, Polyhedron& mesh) {
    Halfedge_handle h = f->halfedge(); int s=0;
    while (h->next()->next()->next() != h) { mesh.split_facet(h, h->next()->next()); if(++s>500) break; }
}

void full_sweep_triangulation(Polyhedron& mesh) {
    std::vector<Facet_handle> facets;
    for (auto f = mesh.facets_begin(); f != mesh.facets_end(); ++f) facets.push_back(f);
    int count = 0, fallbacks = 0;
    for (auto f : facets) {
        if (!triangulate_facet_sweepline(f, mesh)) { triangulate_facet_fan(f, mesh); fallbacks++; }
        if (++count % 100 == 0) std::cout << "." << std::flush;
    }
    facets.clear(); for (auto f = mesh.facets_begin(); f != mesh.facets_end(); ++f) facets.push_back(f);
    for (auto f : facets) if (!f->is_triangle()) triangulate_monotone_cleanup(f, mesh);
}

void save_to_file(const Polyhedron& mesh, const std::string& filename) {
    std::ofstream out(filename); out << mesh;
    std::cout << "Saved to " << filename << std::endl;
}

int main() {
    while (true) {
        Polyhedron mesh;
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "1. Generate New Mesh" << std::endl;
        std::cout << "2. Load OFF File" << std::endl;
        std::cout << "3. Exit" << std::endl;
        std::cout << "Select: ";
        int choice; if (!(std::cin >> choice)) return 0;
        
        if (choice == 1) {
            int n; std::cout << "Points (e.g. 5000): "; std::cin >> n;
            generate_star_mesh(n, mesh);
        } else if (choice == 2) {
            std::string fn; std::cout << "Filename: "; std::cin >> fn;
            std::ifstream in(fn); if (!in || !(in >> mesh)) return 1;
        } else {
            break;
        }
    
        CGAL::draw(mesh); 
        
        full_sweep_triangulation(mesh);
        
        save_to_file(mesh, "result.off");
        CGAL::draw(mesh); 
    }
    return 0;
}
