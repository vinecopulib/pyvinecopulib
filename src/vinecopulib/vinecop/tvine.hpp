#pragma once

#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_optimization.hpp>
#include <wdm/eigen.hpp>

// TODO: weights

namespace vinecopulib {

namespace tools_stl {

template<class T>
std::vector<T> span(std::vector<T> x, size_t start, size_t len)
{
    x.erase(x.begin(), x.begin() + std::min(x.size(), start));
    x.resize(std::min(x.size(), len));
    return x;
}

template<class T>
std::vector<T> rev(std::vector<T> x)
{
    tools_stl::reverse(x);
    return x;
}

template<typename T>
T min_vec(std::vector<T>& x)
{
    auto i = std::min_element(x.begin(), x.end());
    return x[std::distance(x.begin(), i)];
}

template<typename T>
T max_vec(std::vector<T>& x)
{
    auto i = std::max_element(x.begin(), x.end());
    return x[std::distance(x.begin(), i)];
}


}

// ------------------------- T-VINE STRUCTURE ---------------

TriangularArray<size_t> get_actual_struct_array(const RVineStructure& structure)
{
    // structure array is in natural order, must convert to actual order
    auto no_array = structure.get_struct_array();
    auto order = structure.get_order();
    for (size_t i = 0; i < structure.get_dim() - 1; i++) {
        for (auto &s : no_array[i]) {
            s = order[s - 1];
        }
    }

    return no_array;
}

class TVineStructure : public RVineStructure {
public:
    TVineStructure() : RVineStructure() {}

    TVineStructure(size_t cs_dim, size_t p) :
        TVineStructure(RVineStructure(tools_stl::seq_int(1, cs_dim)), p,
                       tools_stl::seq_int(1, cs_dim),
                       tools_stl::seq_int(1, cs_dim))
    {}

    TVineStructure(const RVineStructure &cs_struct,
                   size_t p,
                   std::vector<size_t> in_vertices,
                   std::vector<size_t> out_vertices) :
        p_(p),
        in_vertices_(in_vertices),
        out_vertices_(out_vertices)
    {
        check_in_out_vertices(cs_struct, in_vertices, out_vertices);
        std::cout << cs_struct << std::endl;
        cs_struct_ = reorder_structure(cs_struct, in_vertices);
        std::cout << cs_struct << std::endl;
        order_ = expand_order(cs_struct_.get_order(), p);
        struct_array_ = build_t_vine_array(cs_struct_, p,
                                           in_vertices, out_vertices);

        RVineStructure new_struct;
        try {
            new_struct = RVineStructure(order_, struct_array_);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "out_vertices are not compatible with cross-sectional structure.");
        }
        d_             = new_struct.get_dim();
        trunc_lvl_     = new_struct.get_trunc_lvl();
        struct_array_  = new_struct.get_struct_array();
        min_array_     = new_struct.get_min_array();
        needed_hfunc1_ = new_struct.get_needed_hfunc1();
        needed_hfunc2_ = new_struct.get_needed_hfunc2();
    }

    size_t get_p() const
    {
        return p_;
    }

    size_t get_cs_dim() const
    {
        return cs_struct_.get_dim();
    }

    std::vector<size_t> get_in_vertices() const
    {
        return in_vertices_;
    }

    std::vector<size_t> get_out_vertices() const
    {
        return out_vertices_;
    }

    RVineStructure get_cs_structure() const
    {
        return cs_struct_;
    }

private:

    void check_in_out_vertices(const RVineStructure &cs_struct,
                               std::vector<size_t> in_vertices,
                               std::vector<size_t> out_vertices) const
    {
        auto d = cs_struct.get_dim();
        if (!tools_stl::is_same_set(in_vertices, tools_stl::seq_int(1, d)))
            throw std::runtime_error(
                "in_vertices must contain numbers 1, ..., cs_dim.");
        if (!tools_stl::is_same_set(out_vertices, tools_stl::seq_int(1, d)))
            throw std::runtime_error(
                "out_vertices must contain numbers 1, ..., cs_dim.");
    }

    std::vector<size_t> expand_order(const std::vector<size_t> &order, size_t p) const
    {
        size_t cs_dim = order.size();
        size_t d = cs_dim * (p + 1);
        std::vector<size_t> new_order(d);
        for (size_t i = 0; i < d; i++) {
            new_order[i] = order[i % cs_dim] + ((d - 1 - i) / cs_dim) * cs_dim;
        }

        return new_order;
    }

    std::vector<size_t> sup_diag(const std::vector<size_t>& diag,
                                 const TriangularArray<size_t>& struct_array,
                                 size_t index,
                                 size_t column)  const
    {
        size_t d = diag.size();
        std::vector<size_t> x(d - 1);

        size_t i = 0;
        while (diag[i] != index)  {
            x[i] = diag[i];
            i++;
        }

        size_t pivot_col = i;
        while (i < d - 1) {
            x[i] = struct_array(d - 2 - i, pivot_col);
            i++;
        }

        return tools_stl::span(tools_stl::rev(x), 0, d - 1 - column);
    }

    RVineStructure reorder_structure(const RVineStructure& structure,
                                     std::vector<size_t> in_vertices) const
    {
        using namespace tools_stl;
        size_t d = structure.get_dim();
        if (structure.get_trunc_lvl() < d - 1) {
            throw std::runtime_error("T-vines cannot be truncated.");
        }

        // structure array is in natural order, must convert to actual order
        auto old_struct = get_actual_struct_array(structure);
        std::cout << old_struct << std::endl;
        auto new_struct = old_struct;

        // prepare objects
        auto old_order = structure.get_order();
        auto new_order = rev(in_vertices);

        // loop through all columns
        for (size_t i = 0; i < d - 1; i++) {
            auto new_column = sup_diag(old_order, old_struct, new_order[i], i);
            auto diag_until = span(new_order, 0, i);
            auto diag_after = span(new_order, i, d - i);
            auto max_col = find_position(in_vertices[i], old_order);
            for (auto o : new_column)
                std::cout << o << ", ";
            std::cout << std::endl;
            for (size_t t = 0; t < new_column.size(); t++) {
                // Check whether an element in this column is already contained in
                // the diagonal to the left. If so, we need to find another node
                // that is connected to the diagonal entry of column i. We search
                // for such an edge in the old structure, but only to the the left of
                // the current column or the column, where the element appeared
                // on the diagonal (included).
                if (is_member(new_column[t], diag_until)) {
                    bool found_node = false;
                    for (size_t j = 0; j <= std::min(i, max_col); j++) {
                        if (new_order[i] == old_struct(t, j)) {
                            if (is_member(old_order[j], diag_after)) {
                                new_column[t] = old_order[j];
                                found_node = true;
                            }
                        } else if (new_order[i] == old_order[j]) {
                            if (is_member(old_struct(t, j), diag_after)) {
                                new_column[t] = old_struct(t, j);
                                found_node = true;
                            }
                        }
                        if (found_node) {
                            // The new entry may already be contained in this
                            // column. We need to check the next rows for that
                            // as well.
                            diag_until = cat(new_column[t], diag_until);
                            break;
                        }
                    }
                }
            }
            new_struct[i] = new_column;
            std::cout << new_struct << std::endl;
        }

        // this must always hold beacuse the first in vertex comes last on the diagonal:
        new_struct(0, d - 2) = in_vertices[0];

        RVineStructure new_rvine;
        try {
            new_rvine = RVineStructure(new_order, new_struct);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "in_vertices are not compatible with cross-sectional structure.");
        }

        return RVineStructure(new_order, new_struct);
    }

    TriangularArray<size_t> build_t_vine_array(
        const RVineStructure& cs_struct,
        size_t p,
        std::vector<size_t> in_vertices,
        std::vector<size_t> out_vertices)
        const
    {
        size_t cs_dim = cs_struct.get_dim();
        size_t d = cs_dim * (p + 1);
        auto diag = cs_struct.get_order();

        RVineStructure new_struct = cs_struct;
        if (diag[cs_dim - 1] != in_vertices[0]) {
            new_struct = reorder_structure(new_struct, in_vertices);
        }

        auto struct_array = get_actual_struct_array(cs_struct);

        TriangularArray<size_t> strct(d);
        // copy cross-sectional structure
        for (size_t i = 0; i < cs_dim - 1; i++) {
            for (size_t j = 0; j < cs_dim - 1 - i; j++) {
                strct(i, j) = struct_array(i, j) + cs_dim * p;
            }
        }

        // fill parallelograms
        std::vector<size_t> par = out_vertices;
        for (size_t lag = 1; lag <= p; lag++) {
            for (size_t i = 0; i < cs_dim; i++) {
                for (size_t j = 0; j < cs_dim; j++) {
                    strct(i + cs_dim * lag - j - 1, j) = par[i] + cs_dim * (p - lag);
                }
            }
        }

        // copy to other lags
        for (size_t lag = 1; lag <= p; lag++) {
            for (size_t j = 0; j < cs_dim; j++) {
                for (size_t i = 0; i < d - 1 - (j + cs_dim * lag); i++) {
                    strct(i, j + cs_dim * lag) = strct(i, j) - cs_dim * lag;
                }
            }
        }

        return strct;
    }

private:
    size_t p_;
    std::vector<size_t> in_vertices_;
    std::vector<size_t> out_vertices_;
    RVineStructure cs_struct_;
};


// ------------------------- SELECTOR ------------------------

inline Eigen::MatrixXd spread_lag(const Eigen::MatrixXd& data, size_t cs_dim)
{
    if (data.rows() < 2) {
        throw std::runtime_error("insufficient number of observations");
    }
    if (data.cols() % cs_dim != 0) {
        throw std::runtime_error("number of columns is not a multiple of cs_dim");
    }
    size_t n = data.rows() - 1;
    Eigen::MatrixXd newdata(n, data.cols() + cs_dim);
    newdata << data.topRows(n), data.rightCols(cs_dim).bottomRows(n);
    return newdata;
}

namespace tools_select {

class TVineSelector {
public:
    TVineSelector(const Eigen::MatrixXd &data,
                  std::vector<size_t> in_vertices,
                  std::vector<size_t> out_vertices) :
        cs_dim_(data.cols()),
        lag_(0),
        in_vertices_(in_vertices),
        out_vertices_(out_vertices),
        data_(data)
    {
        check_in_out_vertices();
    }

    TVineSelector(const Eigen::MatrixXd &data) :
        cs_dim_(data.cols()),
        lag_(0),
        data_(data)
    {
    }

    std::vector<size_t> get_in_vertices() const
    {
        return in_vertices_;
    }

    std::vector<size_t> get_out_vertices() const
    {
        return out_vertices_;
    }

    RVineStructure get_cs_structure() const
    {
        return cs_struct_;
    }

protected:
    void duplicate_vertex(size_t v, VineTree& tree)
    {
        auto v_new = boost::add_vertex(tree);
        auto shift = [this] (std::vector<size_t> index) {
            for (auto &i : index)
                i = i + cs_dim_ * lag_;
            return index;
        };

        // copy structure information
        tree[v_new].conditioned = shift(tree[v].conditioned);
        tree[v_new].conditioning = shift(tree[v].conditioning);
        tree[v_new].all_indices = shift(tree[v].all_indices);
        tree[v_new].prev_edge_indices = shift(tree[v].prev_edge_indices);

        // copy data and remove rows
        size_t n = tree[v].hfunc1.rows() - 1;
        tree[v_new].hfunc1 = tree[v].hfunc1.bottomRows(n);
        tree[v].hfunc1.conservativeResize(n);
        if (tree[v].hfunc2.size() > 1) {
            tree[v_new].hfunc2 = tree[v].hfunc2.bottomRows(n);
            tree[v].hfunc2.conservativeResize(n);
        }
    }

    void duplicate_edge(EdgeIterator e, VineTree& tree)
    {
        size_t v1 = boost::source(e, tree);
        size_t v2 = boost::target(e, tree);
        auto e_new = boost::add_edge(v1 + cs_dim_ * lag_, v2 + cs_dim_ * lag_, tree);
        tree[e_new.first].pair_copula = tree[e].pair_copula;
        tree[e_new.first].fit_id = tree[e].fit_id;
        auto shift = [this] (std::vector<size_t> index) {
            for (auto &i : index)
                i = i + cs_dim_ * lag_;
            return index;
        };
        tree[e_new.first].conditioned = shift(tree[e].conditioned);
        tree[e_new.first].conditioning = shift(tree[e].conditioning);
        tree[e_new.first].all_indices = shift(tree[e].all_indices);
    }

    void check_in_out_vertices() const
    {
        auto d = cs_dim_;
        if (!tools_stl::is_same_set(in_vertices_, tools_stl::seq_int(1, d)))
            throw std::runtime_error(
                "in_vertices must contain numbers 1, ..., cs_dim.");
        if (!tools_stl::is_same_set(out_vertices_, tools_stl::seq_int(1, d)))
            throw std::runtime_error(
                "out_vertices must contain numbers 1, ..., cs_dim.");
    }

    void check_controls(const FitControlsVinecop &controls)
    {
        if (controls.get_select_trunc_lvl()) {
            throw std::runtime_error("Cannot select truncation level for T-vines.");
        }
        if (controls.get_trunc_lvl() < std::numeric_limits<int>::max()) {
            throw std::runtime_error("T-vines cannot be truncated.");
        }
    }

    size_t cs_dim_;
    size_t lag_;
    std::vector<size_t> in_vertices_;
    std::vector<size_t> out_vertices_;
    Eigen::MatrixXd data_;
    RVineStructure cs_struct_;
};


class TVineStructureSelector : public StructureSelector, public TVineSelector {
public:
    TVineStructureSelector(const Eigen::MatrixXd &data,
                           const FitControlsVinecop &controls) :
        StructureSelector(data, controls),
        TVineSelector(data)
    {
        check_controls(controls);
    }

    ~TVineStructureSelector() = default;

    void add_lag()
    {
        controls_.set_trunc_lvl(std::numeric_limits<size_t>::max());
        lag_++;
        d_ += cs_dim_;
        data_ = spread_lag(data_, cs_dim_);
        trees_.resize(d_);
        trees_opt_.resize(d_);

        vine_struct_ = RVineStructure(tools_stl::seq_int(1, d_), 1, false);
        trees_[0] = make_base_tree(data_);
        // add vertices and edges for lagged variable
        for (size_t t = 1; t < trees_.size(); t++) {
            auto old_tree = trees_[t];
            auto new_tree = edges_as_vertices(trees_[t - 1]);
            duplicate_edges(old_tree, new_tree, t);
            trees_opt_[t] = new_tree;
            min_spanning_tree(new_tree);
            add_edge_info(new_tree);
            select_pair_copulas(new_tree, trees_opt_[t]);

            trees_[t] = new_tree;
            if (controls_.get_show_trace()) {
                std::cout << "** Tree: " << t - 1 << std::endl;
                print_pair_copulas_of_tree(t - 1);
            }
        }
    }

    void duplicate_edges(VineTree& old_tree, VineTree& new_tree, size_t t)
    {
        // copy existing edges
        for (auto e : boost::edges(old_tree)) {
            size_t v0 = boost::source(e, old_tree);
            size_t v1 = boost::target(e, old_tree);
            auto e_new = boost::add_edge(v0, v1, new_tree).first;
            new_tree[e_new].pair_copula = old_tree[e].pair_copula;
            new_tree[e_new].fit_id = old_tree[e].fit_id;
            new_tree[e_new].conditioned = old_tree[e].conditioned;
            new_tree[e_new].conditioning = old_tree[e].conditioning;
            new_tree[e_new].all_indices = old_tree[e].all_indices;
        }

        if (t >= (lag_ - 1) * cs_dim_)
            add_allowed_connections(new_tree, t);

        // add edges for new lag
        auto shift = [this] (std::vector<size_t> index) {
            for (auto &i : index)
                i = i + cs_dim_ ;
            return index;
        };
        long int ei = 0;
        long int e_start = (lag_ - 1) * cs_dim_ - t;
        long int ei_max = lag_ * cs_dim_ - t;
        for (auto e : boost::edges(new_tree)) {
            if (ei++ < e_start)
                continue;
            if (ei > ei_max)
                continue;

            size_t v0 = boost::source(e, old_tree);
            size_t v1 = boost::target(e, old_tree);
            auto e_new = boost::add_edge(v0 + cs_dim_,
                                         v1 + cs_dim_,
                                         new_tree).first;
            new_tree[e_new].pair_copula = old_tree[e].pair_copula;
            new_tree[e_new].fit_id = old_tree[e].fit_id;
            new_tree[e_new].conditioned = shift(old_tree[e].conditioned);
            new_tree[e_new].conditioning = shift(old_tree[e].conditioning);
            new_tree[e_new].all_indices = shift(old_tree[e].all_indices);
        }
    }

    void add_allowed_connections(VineTree& tree, size_t t)
    {
        auto add_edge = [&tree, this](size_t v0, size_t v1) {
            auto crit = std::fabs(wdm::wdm(get_pc_data(v0, v1, tree), "hoeffd")(0, 1));
            auto w = 1 - crit;
            auto e = boost::add_edge(v0, v1, w, tree).first;
            tree[e].weight = w;
            tree[e].crit = 1 - w;
        };
        auto nv = boost::num_vertices(tree);
        if (t == 1) {
            for (size_t v0 = 0; v0 < cs_dim_ - t + 1; ++v0) {
                for (size_t v1 = cs_dim_ - t + 1; v1 < nv; ++v1) {
                    add_edge(v0, v1);
                }
            }
        } else {
            long int vstart = 0;
            if (t < cs_dim_) {
                vstart = cs_dim_ - t + 1;
            }
            for (size_t v0 = vstart; v0 < std::min(cs_dim_, nv); ++v0) {
                for (size_t v1 = 0; v1 < nv; ++v1) {
                    if (v0 == v1)
                        continue;
                    find_common_neighbor(v0, v1, tree);
                    if (find_common_neighbor(v0, v1, tree) > -1) {
                        if (!boost::edge(v1, v0, tree).second) {
                            add_edge(v0, v1);
                        }
                    }
                }
            }
        }
    }

    void finalize_tvine(size_t trunc_lvl)
    {
        using namespace tools_stl;
        trees_opt_ = trees_;
        pair_copulas_ = make_pair_copula_store(d_, trunc_lvl);
        TriangularArray<size_t> mat(d_, trunc_lvl);
        std::vector<size_t> order(d_);

        if (d_ > cs_dim_) {
            size_t min_c, max_c, min_D, max_D;
            for (size_t t = 0; t < cs_dim_; t++) {
                for (auto e : boost::edges(trees_[t + 1])) {
                    min_c = min_vec(trees_[t + 1][e].conditioned);
                    max_c = max_vec(trees_[t + 1][e].conditioned);
                    if (max_c / cs_dim_ != 1)
                        continue;
                    if (t > 0) {
                        min_D = min_vec(trees_[t + 1][e].conditioning);
                        max_D = max_vec(trees_[t + 1][e].conditioning);
                        if (max_D / cs_dim_ == 0)
                            out_vertices_[t] = min_c;
                        if ((min_D / cs_dim_ == 1) && (min_c  / cs_dim_ == 0))
                            in_vertices_[t] = max_c % cs_dim_;
                    } else {
                        if (min_c  / cs_dim_ == 0) {
                            out_vertices_[t] = min_c;
                            in_vertices_[t] = max_c % cs_dim_;
                        }
                    }
                }
            }
        }

        if (trunc_lvl > 0) {
            std::vector <size_t> ned_set;
            std::vector <size_t> ning_set;

            // fill matrix column by column
            for (size_t col = 0; col < d_ - 1; ++col) {
                tools_interface::check_user_interrupt();
                // matrix above trunc_lvl is left empty
                size_t t =
                    std::max(std::min(trunc_lvl, d_ - 1 - col), static_cast<size_t>(1));
                // start with highest tree in this column
                for (auto e : boost::edges(trees_[t])) {

                    // find an edge that contains a leaf
                    size_t v0 = boost::source(e, trees_[t]);
                    size_t v1 = boost::target(e, trees_[t]);
                    size_t min_deg = std::min(boost::out_degree(v0, trees_[t]),
                                              boost::out_degree(v1, trees_[t]));
                    if (min_deg > 1) {
                        continue;  // not a leaf
                    }

                    // check if edge contains in_vertex
                    ned_set = trees_[t][e].conditioned;
                    if ((in_vertices_[t % cs_dim_] != (ned_set[0] % cs_dim_)) 
                        && (in_vertices_[t % cs_dim_] != (ned_set[1] % cs_dim_)))
                        continue;

                    // find position of in_vertex in the edge
                    size_t pos = (ned_set[1] % cs_dim_ == in_vertices_[t % cs_dim_]);
                    if (pos == 1) {
                        trees_[t][e].pair_copula.flip();
                    }

                    // fill diagonal entry with in vertex
                    order[col] = ned_set[pos];

                    // entry in row t-1 is other index of the edge
                    mat(t - 1, col) = ned_set[1 - pos];

                    // assign fitted pair copula to appropriate entry, see
                    // `Vinecop::get_pair_copula()`.
                    if (trunc_lvl > 0 ) {
                        pair_copulas_[t - 1][col] = trees_[t][e].pair_copula;
                    }

                    // initialize running set with full conditioning set of this edge
                    ning_set = trees_[t][e].conditioning;

                    // remove edge (must not be reused in another column!)
                    boost::remove_edge(v0, v1, trees_[t]);
                    break;
                }

                // fill column bottom to top
                for (size_t k = 1; k < t; ++k) {
                    auto check_set = cat(order[col], ning_set);
                    for (auto e : boost::edges(trees_[t - k])) {
                        // search for an edge in lower tree that shares all
                        // indices in the conditioning set + diagonal entry
                        if (!is_same_set(trees_[t - k][e].all_indices, check_set)) {
                            continue;
                        }
                        // found suitable edge ->
                        // next matrix entry is conditioned variable of new edge
                        // that's not equal to the diagonal entry of this column
                        auto e_new = trees_[t - k][e];
                        ptrdiff_t pos = (order[col] == e_new.conditioned[1]);
                        if (pos == 1) {
                            e_new.pair_copula.flip();
                        }
                        mat(t - k - 1, col) = e_new.conditioned[std::abs(1 - pos)];

                        // assign fitted pair copula to appropriate entry, see
                        // Vinecop::get_pair_copula().
                        pair_copulas_[t - 1 - k][col] = e_new.pair_copula;

                        // start over with conditioned set of next edge
                        ning_set = e_new.conditioning;

                        // remove edge (must not be reused in another column!)
                        size_t v0 = boost::source(e, trees_[t - k]);
                        size_t v1 = boost::target(e, trees_[t - k]);
                        boost::remove_edge(v0, v1, trees_[t - k]);
                        break;
                    }
                }

            }

            // The last column contains a single element which must be different
            // from all other diagonal elements. Based on the properties of an
            // R-vine matrix, this must be the element next to it.
            order[d_ - 1] = mat(0, d_ - 2);

            // change to user-facing format
            // (variable index starting at 1 instead of 0)
            for (size_t i = 0; i < std::min(d_ - 1, trunc_lvl); ++i) {
                for (size_t j = 0; j < d_ - i - 1; ++j) {
                    mat(i, j) += 1;
                }
            }
            for (size_t i = 0; i < d_; i++) {
                order[i] += 1;
            }
            for (size_t i = 0; i < cs_dim_; i++) {
                in_vertices_[i] += 1;
                out_vertices_[i] += 1;
            }
        } else {
            // order doesn't matter for truncated
            order = tools_stl::seq_int(1, d_);
        }

        // return as RVineStructure
        vine_struct_ = RVineStructure(order, mat);
    }

    void finalize(size_t trunc_lvl)
    {
        if (d_ == cs_dim_) {
            StructureSelector::finalize(trunc_lvl);
            trees_ = trees_opt_;
            auto mat = vine_struct_.get_matrix();
            cs_struct_ = RVineStructure(mat.block(0, d_ - cs_dim_, cs_dim_, cs_dim_));
            in_vertices_ = tools_stl::rev(cs_struct_.get_order());
            out_vertices_ = in_vertices_;
        } else {
            finalize_tvine(trunc_lvl);
            auto mat = vine_struct_.get_matrix();
            cs_struct_ = RVineStructure(mat.block(0, d_ - cs_dim_, cs_dim_, cs_dim_));
        }
    }

protected:
    double compute_fit_id(const EdgeProperties& e) override
    {
        size_t min_c = std::min(e.conditioned[0], e.conditioned[1]);
        size_t max_c = std::max(e.conditioned[0], e.conditioned[1]);
        size_t lagdiff = (max_c  - min_c) / cs_dim_;
        return (min_c % cs_dim_) +
                    (max_c % cs_dim_) * std::pow(cs_dim_, 2) +
                    lagdiff * std::pow(cs_dim_, 4);
    }


};

class TVineFamilySelector : public FamilySelector, public TVineSelector {
public:
    TVineFamilySelector(const Eigen::MatrixXd &data,
                        const RVineStructure &cs_struct,
                        const FitControlsVinecop &controls,
                        std::vector<size_t> in_vertices,
                        std::vector<size_t> out_vertices)  :
        FamilySelector(data, cs_struct, controls),
        TVineSelector(data, in_vertices, out_vertices)
    {
        check_controls(controls);
        cs_struct_ = TVineStructure(cs_struct, 0, in_vertices, out_vertices);
    }

    ~TVineFamilySelector() {}

    void select_tree(size_t t) override
    {
        auto new_tree = edges_as_vertices(trees_[t]);
        remove_edge_data(trees_[t]);
        add_allowed_edges(new_tree);
        if (boost::num_vertices(new_tree) > 0) {
            add_edge_info(new_tree);       // for pc estimation and next tree
            if (controls_.get_selection_criterion() == "mbicv") {
                // adjust prior probability to tree level
                controls_.set_psi0(std::pow(psi0_, t + 1));
            }
            if (trees_opt_.size() > t + 1) {
                select_pair_copulas(new_tree, trees_opt_[t + 1]);
            } else {
                select_pair_copulas(new_tree);
            }
        }
        // make sure there is space for new tree
        trees_.resize(t + 2);
        trees_[t + 1] = new_tree;
    }

    void add_lag()
    {
        lag_++;
        d_ += cs_dim_;

        // add vertices and edges for lagged variable
        for (size_t t = 1; t < trees_.size(); t++) {
            auto old_tree = trees_[t];
            for (auto v : boost::vertices(old_tree))
                duplicate_vertex(v, trees_[t]);
            for (auto e : boost::edges(old_tree))
                duplicate_edge(e, trees_[t]);
        }

        // update trees and structure
        trees_opt_ = trees_;
        trees_ = std::vector<VineTree>(1);
        vine_struct_ = TVineStructure(cs_struct_, lag_, in_vertices_, out_vertices_);
        data_ = spread_lag(data_, cs_dim_);
        controls_.set_trunc_lvl(std::numeric_limits<size_t>::max());
    }

    Eigen::MatrixXd data()
    {
        return data_;
    }

protected:

    void finalize(size_t)
    {
        pair_copulas_ = make_pair_copula_store(d_, d_);
        std::vector<size_t> order = vine_struct_.get_order();
        size_t d = order.size();
        auto struct_array = get_actual_struct_array(vine_struct_);
        for (size_t t = 1; t < trees_.size(); t++) {
            for (auto e : boost::edges(trees_[t])) {
                std::vector<size_t> check_set(2);
                for (size_t k = 0; k < d; k++) {
                    // conditioned starts at 0 -> substract -1
                    check_set = {order[k] - 1, struct_array(t - 1, k) - 1};
                    if (tools_stl::is_same_set(trees_[t][e].conditioned, check_set)) {
                        pair_copulas_[t - 1][k] = trees_[t][e].pair_copula;
                        break;
                    }
                }
            }
        }
    }

    double compute_fit_id(const EdgeProperties& e) override
    {
        size_t min_c = std::min(e.conditioned[0], e.conditioned[1]);
        size_t max_c = std::max(e.conditioned[0], e.conditioned[1]);
        return (min_c % cs_dim_) * cs_dim_ * 10 + (max_c % cs_dim_);
    }

    void flip_pcs(std::vector<VineTree>& trees)
    {
        std::vector<size_t> order = cs_struct_.get_order();
        size_t d = order.size();
        auto struct_array = get_actual_struct_array(cs_struct_);
        for (size_t t = 1; t < trees.size(); t++) {
            for (auto e : boost::edges(trees[t])) {
                std::vector<size_t> check_set(2);
                for (size_t k = 0; k < d; k++) {
                    // conditioned starts at 0 -> substract -1
                    check_set = {order[k] - 1, struct_array(t - 1, k) - 1};
                    if (tools_stl::is_same_set(trees_[t][e].conditioned, check_set)) {
                        break;
                    }
                }
                if (trees[t][e].conditioned[0] != check_set[0]) {
                    trees[t][e].pair_copula.flip();
                }
            }
        }
    }

};

} // end tools_select


// --------------------- T-VINE ----------------------------------




class TVine : public Vinecop {
public:
    TVine(size_t cs_dim, size_t p) :
        TVine(tools_stl::seq_int(1, cs_dim), p,
              tools_stl::rev(tools_stl::seq_int(1, cs_dim)),
              tools_stl::rev(tools_stl::seq_int(1, cs_dim)))
    {}

    TVine(const RVineStructure &cs_struct,
          size_t p,
          std::vector<size_t> in_vertices,
          std::vector<size_t> out_vertices) :
        cs_dim_(cs_struct.get_dim()),
        p_(p),
        in_vertices_(in_vertices),
        out_vertices_(out_vertices),
        tvine_struct_(TVineStructure(cs_struct, p, in_vertices, out_vertices))
    {
        threshold_ = 0.0;
        loglik_ = NAN;
	    d_ = tvine_struct_.get_dim();
        vine_struct_ = tvine_struct_;
        pair_copulas_ = make_pair_copula_store(d_);
    }

    TVine(const std::vector<std::vector<Bicop>> &pair_copulas,
          const RVineStructure &cs_struct,
          size_t p,
          std::vector<size_t> in_vertices,
          std::vector<size_t> out_vertices) :
         TVine(cs_struct, p, in_vertices, out_vertices)
    {
        pair_copulas_ = pair_copulas;
    }

    size_t get_p() const
    {
        return p_;
    }

    size_t get_cs_dim() const
    {
        return cs_dim_;
    }

    std::vector<size_t> get_in_vertices() const
    {
        return in_vertices_;
    }

    std::vector<size_t> get_out_vertices() const
    {
        return out_vertices_;
    }

    RVineStructure get_cs_structure() const
    {
        return tvine_struct_.get_cs_structure();
    }

    TVineStructure get_tvine_structure() const
    {
        return tvine_struct_;
    }

    void select_families(const Eigen::MatrixXd &data,
                         const FitControlsVinecop &controls = FitControlsVinecop())
    {
        tools_eigen::check_if_in_unit_cube(data);
        check_data_dim(data);

        if (vine_struct_.get_trunc_lvl() > 0) {
            tools_select::TVineFamilySelector selector(
                data,
                tvine_struct_.get_cs_structure(),
                controls,
                in_vertices_,
                out_vertices_);

            selector.select_all_trees(data);
            for (size_t lag = 1; lag <= p_; lag++) {
                selector.add_lag();
                selector.select_all_trees(selector.data());
            }

            finalize_fit(selector);
            loglik_ = this->loglik(data);
        }
    }

    void select_all(const Eigen::MatrixXd &data,
                    const FitControlsVinecop &controls = FitControlsVinecop())
    {
        tools_eigen::check_if_in_unit_cube(data);
        check_data_dim(data);

        tools_select::TVineStructureSelector selector(data, controls);
        selector.select_all_trees(data);
        for (size_t lag = 1; lag <= p_; lag++) {
            selector.add_lag();
        }
        finalize_fit(selector);
        loglik_ = this->loglik(data);
    }

    Eigen::MatrixXd simulate(const size_t n,
                             const bool qrng = false,
                             const size_t num_threads = 1,
                             const std::vector<int>& seeds = std::vector<int>()) const
    {
        Eigen::MatrixXd U = my_simulate_uniform(n, cs_dim_, qrng, seeds);

        // initialize first p + 1 lags
        Eigen::MatrixXd sim(n, cs_dim_);
        Eigen::MatrixXd Ui(1, d_);
        for (size_t i = 0; i <= p_; i++) {
            Ui.row(0).segment(i * cs_dim_, cs_dim_) = U.row(i);
        }
        Eigen::MatrixXd V = inverse_rosenblatt(Ui, num_threads);
        for (size_t i = 0; i <= p_; i++) {
            sim.row(i) = V.block(0, i * cs_dim_, 1, cs_dim_);
        }

        // simulate conditional on previous observations
        for (size_t i = p_ + 1; i < n; i++) {
            Ui.leftCols(d_ - cs_dim_).swap(Ui.rightCols(d_ - cs_dim_));
            Ui.rightCols(cs_dim_) = U.row(i);
            sim.row(i) = inverse_rosenblatt(Ui, num_threads).rightCols(cs_dim_);
        }

        return sim;
    }

    Eigen::MatrixXd simulate_conditional(size_t n,
        const Eigen::MatrixXd& data,
        const bool qrng = false,
        const size_t num_threads = 1,
        const std::vector<int>& seeds = std::vector<int>())
    {
        check_cond_data(data);

        Eigen::MatrixXd U(n, d_);
        if (p_ > 0) {
            U.leftCols(d_ - cs_dim_) = get_last_cpits(data).replicate(n, 1);
        }
        U.rightCols(cs_dim_) = my_simulate_uniform(n, cs_dim_, qrng, seeds);

        return inverse_rosenblatt(U, num_threads).rightCols(cs_dim_);
    }

    Eigen::MatrixXd simulate_ahead(size_t n_ahead,
                                   const Eigen::MatrixXd& data,
                                   const bool qrng = false,
                                   const std::vector<int>& seeds = std::vector<int>())
    {
        check_cond_data(data);

        Eigen::MatrixXd U(n_ahead + p_, cs_dim_);
        U.bottomRows(n_ahead) = my_simulate_uniform(n_ahead, cs_dim_, qrng, seeds);
        Eigen::MatrixXd Ui(1, d_);

        // initialize first p lags
        if (p_ > 0) {
            Ui.rightCols(d_ - cs_dim_) = get_last_cpits(data);
        }

        // simulate conditional on previous observations
        Eigen::MatrixXd sim(n_ahead, cs_dim_);
        for (size_t i = 0; i < n_ahead; i++) {
            Ui.leftCols(d_ - cs_dim_).swap(Ui.rightCols(d_ - cs_dim_));
            Ui.rightCols(cs_dim_) = U.row(i + p_);
            sim.row(i) = inverse_rosenblatt(Ui).rightCols(cs_dim_);
        }

        return sim;
    }

    Eigen::VectorXd pdf(const Eigen::MatrixXd &, const size_t = 1) const
    {
        throw std::runtime_error("pdf not meaningful for T-vines; use loglik().");
    }

    double loglik(const Eigen::MatrixXd &u, const size_t num_threads = 1)
    {
        if (static_cast<size_t>(u.cols()) != cs_dim_)
            throw std::runtime_error("dimension of data and model don't match.");
        size_t n = u.rows();

        // iid model, can return loglik directly
        if ((p_ == 0) || (n == 1)) {
            vine_struct_ = tvine_struct_.get_cs_structure();
            return Vinecop::loglik(u, num_threads);
        }

        // first compute substraction component (otherwise some contributions
        // are counted twice)
        size_t p_tmp = std::min(n, p_) - 1;
        vine_struct_ = TVineStructure(tvine_struct_.get_cs_structure(),
                                      p_tmp, in_vertices_, out_vertices_);
        d_ = cs_dim_ * (1 + p_tmp);
        auto u_spr = u;
        for (size_t lag = 0; lag < p_tmp; ++lag) {
            u_spr = spread_lag(u_spr, cs_dim_);
        }

        n = u_spr.rows();
        double ll = 0.0;
        if (n > 2) {
            ll -= Vinecop::loglik(u_spr.bottomRows(n - 1).topRows(n - 2),
                                        num_threads);
        } else {
            ll -= Vinecop::loglik(u_spr.bottomRows(n - 1), num_threads);
        }

        // add loglik *as if* it was iid with cs_dim * (1 + p) vars
        u_spr = spread_lag(u_spr, cs_dim_);
        vine_struct_ = TVineStructure(tvine_struct_.get_cs_structure(),
                                      p_, in_vertices_, out_vertices_);
        d_ = cs_dim_ * (1 + p_);
        ll += Vinecop::loglik(u_spr, num_threads);

        return ll;
    }

    Eigen::VectorXd cond_cdf(const Eigen::MatrixXd &u, size_t conditioned,
                             const size_t num_threads = 1) const
    {
        if (static_cast<size_t>(u.cols()) != cs_dim_)
            throw std::runtime_error("dimension of data and model don't match.");
        if (static_cast<size_t>(u.rows()) <= p_)
            throw std::runtime_error("insufficient number of time points.");

        auto v = u;
        for (size_t lag = 0; lag < p_; ++lag) {
            v = spread_lag(v, cs_dim_);
            conditioned += cs_dim_;
        }

        size_t n = v.rows();
        Eigen::VectorXd seq = Eigen::VectorXd::LinSpaced(100, 1e-10, 1 - 1e-10);
        Eigen::MatrixXd vv;
        Eigen::VectorXd out(n);
        for (size_t i = 0; i < n; ++i) {
            vv = v.row(i).replicate(100, 1);
            vv.col(conditioned) = seq;
            auto pdf = Vinecop::pdf(vv, num_threads);
            out(i) = pdf.head(std::ceil(v(i, conditioned) * 100)).sum();
            out(i) /= pdf.sum();
        }
        return out;
    }

    Eigen::VectorXi get_num_pars()
    {
        Eigen::VectorXi nums(cs_dim_ * cs_dim_ * p_ + cs_dim_ * (cs_dim_ - 1) / 2);
        size_t i = 0;
        for (size_t t = 0; t < cs_dim_ * (1 + p_) - 1; ++t) {
            for (size_t e = 0; e < cs_dim_; ++e) {
                if (e < pair_copulas_[t].size()) {
                    if (pair_copulas_[t][e].get_family() == BicopFamily::tll) {
                        nums(i++) = 0;
                    } else {
                        nums(i++) = pair_copulas_[t][e].get_parameters().size();
                    }
                }
            }
        }

        return nums;
    }

    Eigen::VectorXd get_pars()
    {
        auto nums = this->get_num_pars();
        Eigen::VectorXd pars(nums.array().sum());
        size_t i = 0, j = 0;
        for (size_t t = 0; t < cs_dim_ * (1 + p_) - 1; ++t) {
            for (size_t e = 0; e < cs_dim_; ++e) {
                if (e < pair_copulas_[t].size()) {
                    if (nums(j) > 0) {
                        pars.segment(i, nums(j)) = Eigen::Map<Eigen::VectorXd>(
                            pair_copulas_[t][e].get_parameters().data(),
                            nums(j)
                        );
                    }
                    i += nums(j++);
                }
            }
        }

        return pars;
    }

    Eigen::VectorXd get_lbs()
    {
        auto nums = this->get_num_pars();
        Eigen::VectorXd pars(nums.array().sum());
        size_t i = 0, j = 0;
        for (size_t t = 0; t < cs_dim_ * (1 + p_) - 1; ++t) {
            for (size_t e = 0; e < cs_dim_; ++e) {
                if (e < pair_copulas_[t].size()) {
                    if (nums(j) > 0) {
                        pars.segment(i, nums(j)) = Eigen::Map<Eigen::VectorXd>(
                            pair_copulas_[t][e].get_parameters_lower_bounds().data(),
                            nums(j)
                        );
                    }
                    i += nums(j++);
                }
            }
        }

        return pars;
    }

    Eigen::VectorXd get_ubs()
    {
        auto nums = this->get_num_pars();
        Eigen::VectorXd pars(nums.array().sum());
        size_t i = 0, j = 0;
        for (size_t t = 0; t < cs_dim_ * (1 + p_) - 1; ++t) {
            for (size_t e = 0; e < cs_dim_; ++e) {
                if (e < pair_copulas_[t].size()) {
                    if (pair_copulas_[t][e].get_family() == BicopFamily::tll) {
                        j++;
                        continue;
                    }
                    if (nums(j) > 0) {
                        pars.segment(i, nums(j)) = Eigen::Map<Eigen::VectorXd>(
                            pair_copulas_[t][e].get_parameters_upper_bounds().data(),
                            nums(j)
                        );
                    }
                    i += nums(j++);
                }
            }
        }

        return pars;
    }


    void set_pars(const Eigen::VectorXd& pars)
    {
        auto nums = this->get_num_pars();
        size_t i = 0, j = 0;
        for (size_t t = 0; t < cs_dim_ * (1 + p_) - 1; ++t) {
            for (size_t e = 0; e < cs_dim_; ++e) {
                if (e >= pair_copulas_[t].size())
                    continue;
                for (size_t lag = 0; lag < p_ + 1; ++lag) {
                    size_t el = e + lag * cs_dim_;
                    if (el < pair_copulas_[t].size()) {
                        pair_copulas_[t][el].set_parameters(pars.segment(i, nums(j)));
                    }
                }
                i += nums(j++);
            }
        }
    }


    typedef struct
    {
        const Eigen::MatrixXd &U; //!< the data.
        vinecopulib::TVine *tvine; //!< a pointer to the vine copula to optimize.
        Eigen::VectorXd pars;  //!< main dependence parameters.
        size_t objective_calls; //!< number of evaluations of the objective.
        const Eigen::VectorXd &weights; //!< weights for the observations.
        double objective_min;  //!< final value of the objective function
    } JointMLEData;


    void fit_joint(const Eigen::MatrixXd& data,
        const FitControlsVinecop &controls = FitControlsVinecop())
    {
        using namespace tools_optimization;

        // Set bounds and starting values
        auto pars = get_pars();
        auto npars = pars.size();
        if (npars == 0) return;

        auto lb = get_lbs();
        auto ub = get_ubs();

        std::function<double(const Eigen::VectorXd&)> objective =
            [this, &data, &controls] (const Eigen::VectorXd& pars) {
            this->set_pars(pars);
            // TODO: weights
            return this->loglik(data, controls.get_num_threads());
        };

        // create optimizer
        Optimizer optimizer;

        // optimize and store result
        pars = optimizer.optimize(pars, lb, ub, objective);

        // update model
        loglik_ = optimizer.get_objective_max();
        set_pars(pars);
    }

protected:
    Eigen::MatrixXd get_last_cpits(const Eigen::MatrixXd& data)
    {
        auto cpits = Eigen::MatrixXd();

        if (p_ > 0) {
            // only most recent observations are used
            Eigen::MatrixXd cond_vals = data.bottomRows(p_);

            // spread past observations into one row with d_ - cs_dim_ columns
            for (size_t lag = 1; lag < p_; lag++) {
                cond_vals = spread_lag(cond_vals, cs_dim_);
            }

            // construct sub-model for last p_ lags
            d_ -= cs_dim_;
            vine_struct_ = TVineStructure(tvine_struct_.get_cs_structure(),
                                          p_ - 1, in_vertices_, out_vertices_);

            // initialize Ui with rosenblatt of past observations
            cpits = rosenblatt(cond_vals);

            // restore original model
            vine_struct_ = tvine_struct_;
            d_ += cs_dim_;
        }

        return cpits;
    }

    Eigen::MatrixXd my_simulate_uniform(
        size_t n, size_t d,
        const bool qrng = false,
        const std::vector<int>& seeds = std::vector<int>()) const
    {
        Eigen::MatrixXd U;
        if (qrng) {
            if (cs_dim_ > 300) {
                U = tools_stats::sobol(n, d, seeds);
            } else {
                U = tools_stats::ghalton(n, d, seeds);
            }
        } else {
            U = tools_stats::simulate_uniform(n, d, seeds);
        }

        return U;
    }

    void finalize_fit(const tools_select::TVineFamilySelector &selector)
    {
        in_vertices_ = selector.get_in_vertices();
        out_vertices_ = selector.get_out_vertices();
        tvine_struct_ = TVineStructure(selector.get_cs_structure(),
                                       p_, in_vertices_, out_vertices_);
        vine_struct_ = tvine_struct_;
        Vinecop::finalize_fit(selector);
    }

    void finalize_fit(tools_select::TVineStructureSelector &selector)
    {
        selector.finalize(std::numeric_limits<size_t>::max());
        in_vertices_ = selector.get_in_vertices();
        out_vertices_ = selector.get_out_vertices();
        Vinecop::finalize_fit(selector);
        tvine_struct_ = TVineStructure(selector.get_cs_structure(),
                                       p_, in_vertices_, out_vertices_);
        vine_struct_ = tvine_struct_;
    }

    void check_data_dim(const Eigen::MatrixXd &data) const
    {
        if (cs_dim_ != static_cast<size_t>(data.cols())) {
            std::stringstream msg;
            msg << "wrong number of columns." << std::endl <<
                "expected: " << cs_dim_ << std::endl <<
                "provided: " << data.cols() << std::endl;
            throw std::runtime_error(msg.str());
        }
    }

    void check_cond_data(const Eigen::MatrixXd &data) const
    {
        check_data_dim(data);
        if (static_cast<size_t>(data.rows()) < p_) {
            std::stringstream msg;
            msg <<
                "need at least p observations to condition on;" << std::endl <<
                "expected: >= " << p_ << std::endl <<
                "actual: " << data.rows() << std::endl;
            throw std::runtime_error(msg.str());
        }
    }


    size_t cs_dim_;
    size_t p_;
    std::vector<size_t> in_vertices_;
    std::vector<size_t> out_vertices_;
    TVineStructure tvine_struct_;
};

}
