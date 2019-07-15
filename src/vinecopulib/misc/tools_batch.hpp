// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vector>

namespace vinecopulib {

namespace tools_batch {

struct Batch {
    size_t begin;
    size_t size;
};

inline std::vector<Batch> create_batches(size_t num_tasks, size_t num_threads)
{
    num_threads = std::max(static_cast<size_t>(1), num_threads);
    std::vector<Batch> batches(std::min(num_tasks, num_threads));
    size_t min_size = num_tasks / num_threads;
    ptrdiff_t rem_size = num_tasks % num_threads;
    for (size_t i = 0, k = 0; i < num_tasks; k++) {
        batches[k] = Batch{i, min_size + (rem_size-- > 0)};
        i += batches[k].size;
    }

    return batches;
}

}

}
