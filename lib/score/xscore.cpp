#include "xscore.hpp"

#include "helper/help.hpp"

inline double get_xs_radius(const Molib::Atom& xs) {
    switch (xs.element().number()) {
        case Molib::Element::C:
            return 1.9;
            break;
        case Molib::Element::N:
            return 1.8;
            break;
        case Molib::Element::O:
            return 1.7;
            break;
        case Molib::Element::S:
            return 2.0;
            break;
        case Molib::Element::P:
            return 2.1;
            break;
        case Molib::Element::F:
            return 1.5;
            break;
        case Molib::Element::Cl:
            return 1.8;
            break;
        case Molib::Element::Br:
            return 2.0;
            break;
        case Molib::Element::I:
            return 2.2;
            break;
        default:
            return 1.2;
            break;
    }
    return 1.2;
}

inline bool xs_is_hydrophobic(const Molib::Atom& xs) {

    switch (xs.element().number()) {

    case Molib::Element::C:

    if (xs.idatm_type() == 17 ||
        xs.idatm_type() == 21)// ||
        //xs.idatm_type() == 22)
        return false;

    for (auto neigh : xs)
        if (neigh.element() == Molib::Element::O ||
            neigh.element() == Molib::Element::N)
            return false;

    case Molib::Element::F:
    case Molib::Element::Cl:
    case Molib::Element::Br:
    case Molib::Element::I:

    return true;

    break;

    default:
    return false;
    break;

    }

    return false;
}

inline bool xs_is_acceptor(const Molib::Atom& xs) {
    if (xs.element() != Molib::Element::N &&
        xs.element() != Molib::Element::O)
        return false;

    auto as = xs.idatm_type();
    if (as == 71 || // IDATM N1+
        as == 73 || // IDATM N2+
        as == 75 || // IDATM N3+
        as == 80 || // IDATM Ng+
        as == 89 || // IDATM O1+
        as == 95 ) {// IDATM Oar+
        return false;
    }

    return true;
}

inline bool xs_is_donor(const Molib::Atom& xs) {

    // Metals
    if (xs.element() == Molib::Element::Mg ||
        xs.element() == Molib::Element::Mn ||
        xs.element() == Molib::Element::Fe ||
        xs.element() == Molib::Element::Zn ||
        xs.element() == Molib::Element::Ca )

        return true;

    if (xs.element() != Molib::Element::N &&
        xs.element() != Molib::Element::O)
        return false;

    int num_h = help::get_info_map(xs.idatm_type_unmask()).substituents - xs.size();

    if (num_h > 0) return true;
    return false;
}

inline double gaussian(double x, double width) {
    const double expo = x/width;
	return std::exp(-(expo*expo));
}

inline double optimal_distance(const Molib::Atom& xs_t1, const Molib::Atom& xs_t2) {
	return get_xs_radius(xs_t1) + get_xs_radius(xs_t2);
}

inline bool xs_donor_acceptor(const Molib::Atom& t1, const Molib::Atom& t2) {
	return xs_is_donor(t1) && xs_is_acceptor(t2);
}

inline bool xs_h_bond_possible(const Molib::Atom& t1, const Molib::Atom& t2) {
	return xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1);
}

inline double slope_step(double x_bad, double x_good, double x) {
	if(x_bad < x_good) {
		if(x <= x_bad) return 0;
		if(x >= x_good) return 1;
	}
	else {
		if(x >= x_bad) return 0;
		if(x <= x_good) return 1;
	}
	return (x - x_bad) / (x_good - x_bad);
}

inline double hydrophobic_eval(const Molib::Atom& t1, const Molib::Atom& t2, double r, double bad, double good) {
	if(xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
		return slope_step(bad, good, r - optimal_distance(t1, t2));
	else return 0;
}

inline double gauss_eval(const Molib::Atom& t1, const Molib::Atom& t2, double r, double offset, double width) {
	return gaussian(r - (optimal_distance(t1, t2) + offset), width);
}

inline double repulsion_eval(const Molib::Atom& t1, const Molib::Atom& t2, double r, double offset) {
	double d = r - (optimal_distance(t1, t2) + offset);
	if(d > 0) 
		return 0;
	return d*d;
}

inline double hydrogen_eval(const Molib::Atom& t1, const Molib::Atom& t2, double r, double bad, double good) {
	if(xs_h_bond_possible(t1, t2))
		return slope_step(bad, good, r - optimal_distance(t1, t2));
	return 0;
}

double Score::vina_xscore(const Molib::Atom::Grid &gridrec, const Molib::Atom::Vec &atoms) {
    double gauss_1 = 0.0;
    double gauss_2 = 0.0;
    double rep = 0.0;
    double hydrophobic = 0.0;
    double hydrogen = 0.0;

    for (size_t i = 0; i < atoms.size(); ++i) {

        const Molib::Atom &atom2 = *atoms[i];

        for (auto &atom1p : gridrec.get_neighbors(atom2.crd(), 8.0)) {
            auto atom1 = *atom1p;
            const double dist = atom1.crd().distance(atom2.crd());

            gauss_1 += gauss_eval( atom1, atom2, dist, 0.0, 0.5 );
            gauss_2 += gauss_eval( atom1, atom2, dist, 3.0, 2.0 );
            rep += repulsion_eval( atom1, atom2, dist, 0);
            hydrophobic += hydrophobic_eval( atom1, atom2, dist, 1.5, 0.5);
            hydrogen += hydrogen_eval(atom1, atom2, dist, 0, -0.7);
        }
    }

    cout << gauss_1 << " " << gauss_2 << " " << rep << " " << hydrophobic << " " << hydrogen << endl;

    return gauss_1 * -0.035579 +
           gauss_2 * -0.005156 +
           rep     *  0.840245 +
           hydrophobic * -0.035069 +
           hydrogen * -0.587439;
}
