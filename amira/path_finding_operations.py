from suffix_tree import Tree


def is_sublist(long_list, sub_list):
    """Check if list is a sublist of long_list."""
    assert isinstance(long_list, list) and isinstance(sub_list, list)
    len_sub = len(sub_list)
    return any(sub_list == long_list[i : i + len_sub] for i in range(len(long_list) - len_sub + 1))


def cluster_adjacent_paths(adjacent_paths):
    # sort the subpaths from longest to shortest
    sorted_paths = sorted([k for k in adjacent_paths], key=len, reverse=True)
    # cluster the sorted paths
    clustered_sub_paths = {}
    for p in sorted_paths:
        list_p = list(p)
        paths_supported = []
        for c in clustered_sub_paths:
            list_c = list(c)
            if is_sublist(list_c, list_p):
                paths_supported.append(c)
        if len(paths_supported) == 0:
            clustered_sub_paths[p] = {p}
        if len(paths_supported) == 1:
            clustered_sub_paths[paths_supported[0]].add(p)
    # choose the shortest subpath in a cluster as the representative
    final_clusters = {}
    for c in clustered_sub_paths:
        final_clusters[min(list(clustered_sub_paths[c]), key=len)] = {
            "longest": max(list(clustered_sub_paths[c]), key=len),
            "all": list(clustered_sub_paths[c]),
        }
    return final_clusters


def find_sublist_indices(main_list, sublist):
    indices = []
    sublist_length = len(sublist)
    # Loop through the main list from start
    for i in range(len(main_list) - sublist_length + 1):
        # Check if the slice of main_list starting at i matches the sublist
        if main_list[i : i + sublist_length] == sublist:
            # Append the start and stop indices to the results
            indices.append((i, i + sublist_length - 1))
    return indices


def construct_suffix_tree(read_nodes):
    rc_reads = {}
    for r in read_nodes:
        if len(set(read_nodes[r])) != 1:
            rc_reads[r + "_reverse"] = list(reversed(read_nodes[r]))
    read_nodes.update(rc_reads)
    return Tree(read_nodes)


def get_suffixes_from_initial_tree(tree, a1):
    suffixes = {}
    for read_id, path in tree.find_all([a1]):
        string_path = str(path).split(" ")
        if string_path[0] == "None":
            continue
        path_list = [int(n) for n in string_path if n != "$"]
        if read_id not in suffixes:
            suffixes[read_id] = path_list
        else:
            if len(path_list) > len(suffixes[read_id]):
                suffixes[read_id] = path_list
    return suffixes


def get_blocks_from_subtree(sub_tree, a2, nodeAnchors):
    block_reads = {}
    block_duplicates = {}
    for read_id, path in sub_tree.find_all([a2]):
        string_path = str(path).split(" ")
        if string_path[0] == "None":
            continue
        path_list = [int(n) for n in string_path if n != "$"]
        assert path_list[0] in nodeAnchors and path_list[-1] in nodeAnchors
        canonical = sorted([path_list, list(reversed(path_list))])[0]
        canonical_tuple = tuple(canonical)
        block_duplicates[canonical_tuple] = False
        if "_reverse" not in read_id:
            if read_id not in block_reads:
                block_reads[read_id] = list(reversed(path_list))
            else:
                if len(path_list) > len(list(reversed(path_list))):
                    block_reads[read_id] = list(reversed(path_list))
    return block_reads, block_duplicates


def get_all_context_options(nodes_on_reads, start, end):
    up_options = set()
    up = nodes_on_reads[:start]
    for i in range(1, len(up) + 1):
        sub_up = tuple(up[-i:])  # Simplified indexing
        up_options.add(sub_up)
    down = nodes_on_reads[end + 1 :]
    down_options = set()
    # Process downstream reads
    for i in range(1, len(down) + 1):
        sub_down = tuple(down[:i])  # Simplified indexing
        down_options.add(sub_down)
    up_options.add(())
    down_options.add(())
    return up_options, down_options


def get_full_path_contexts(positions_of_path, contexts, reads, read_id, block_reads):
    start, end = positions_of_path[0]
    up_options, down_options = get_all_context_options(reads[read_id], start, end)
    canonical = sorted([block_reads[read_id], list(reversed(block_reads[read_id]))])[0]
    canonical_tuple = tuple(canonical)
    if canonical == block_reads[read_id]:
        if canonical_tuple not in contexts:
            contexts[canonical_tuple] = {"upstream": set(), "downstream": set()}
        contexts[canonical_tuple]["upstream"].update(up_options)
        contexts[canonical_tuple]["downstream"].update(down_options)
    else:
        if canonical_tuple not in contexts:
            contexts[canonical_tuple] = {"upstream": set(), "downstream": set()}
        rv_up = set()
        for u in up_options:
            rv_up.add(tuple(reversed(list(u))))
        rv_down = set()
        for d in down_options:
            rv_down.add(tuple(reversed(list(d))))
        contexts[canonical_tuple] = {"upstream": rv_down, "downstream": rv_up}


def orient_nodes_on_read(nodes_on_read, start, end, AMR_indices):
    if AMR_indices[start + 1] == 1:
        nodes_on_read = nodes_on_read
        new_start = start
        new_end = end
    elif AMR_indices[start - 1] == 1:
        nodes_on_read = list(reversed(nodes_on_read))
        new_start = len(nodes_on_read) - 1 - end
        new_end = len(nodes_on_read) - 1 - start
    else:
        nodes_on_read = list(reversed(nodes_on_read))
        new_start = len(nodes_on_read) - 1 - end
        new_end = len(nodes_on_read) - 1 - start
    return new_start, new_end, nodes_on_read


def get_start_stop_indices(binary_list):
    indices = []
    start = None
    for i, value in enumerate(binary_list):
        if value == 1 and start is None:
            start = i  # Start of a sequence
        elif value == 0 and start is not None:
            indices.append((start, i - 1))
            start = None
    # Handle case where the sequence of 1s goes to the end of the list
    if start is not None:
        indices.append((start, len(binary_list) - 1))
    return indices


def get_full_paths(tree, reads, nodeAnchors, threshold):
    # Store full and partial blocks using the suffix tree
    full_blocks = {}
    # Iterate through the anchors
    for a1 in nodeAnchors:
        # Extract paths and build a subtree
        suffixes = get_suffixes_from_initial_tree(tree, a1)
        reversed_suffixes = {}
        for read in suffixes:
            reversed_suffixes[read] = list(reversed(suffixes[read]))
        sub_tree = Tree(reversed_suffixes)
        process_anchors(sub_tree, nodeAnchors, a1, full_blocks, reads, tree, threshold)
    # Filter and return the blocks
    filtered_blocks = filter_blocks(full_blocks)
    return filtered_blocks


def process_anchors(sub_tree, nodeAnchors, a1, full_blocks, reads, tree, threshold):
    for a2 in nodeAnchors:
        if a1 != a2:
            block_reads, block_duplicates = get_blocks_from_subtree(sub_tree, a2, nodeAnchors)
            contexts = generate_contexts(block_reads, block_duplicates, reads)
            generate_full_paths(contexts, block_duplicates, full_blocks, tree, threshold)


def generate_contexts(block_reads, block_duplicates, reads):
    contexts = {}
    for read_id in block_reads:
        positions_of_path = find_sublist_indices(reads[read_id], block_reads[read_id])
        assert len(positions_of_path) > 0
        canonical = get_canonical_representation(block_reads[read_id])
        canonical_tuple = tuple(canonical)
        block_duplicates = update_duplicates(block_duplicates, canonical_tuple, positions_of_path)
        if len(positions_of_path) == 1:
            get_full_path_contexts(positions_of_path, contexts, reads, read_id, block_reads)
    return contexts


def get_canonical_representation(block_read):
    return sorted([block_read, list(reversed(block_read))])[0]


def update_duplicates(block_duplicates, canonical_tuple, positions_of_path):
    if len(positions_of_path) > 1:
        block_duplicates[canonical_tuple] = True
    return block_duplicates


def generate_full_paths(contexts, block_duplicates, full_blocks, tree, threshold):
    for c in contexts:
        if block_duplicates[c] is False:
            upstream_clusters = cluster_adjacent_paths(contexts[c]["upstream"])
            downstream_clusters = cluster_adjacent_paths(contexts[c]["downstream"])
            full_paths = build_full_paths(upstream_clusters, downstream_clusters, c)
            update_full_blocks(full_paths, tree, threshold, full_blocks, c)


def build_full_paths(upstream_clusters, downstream_clusters, c):
    return [u + c + d for u in upstream_clusters for d in downstream_clusters]


def update_full_blocks(full_paths, tree, threshold, full_blocks, c):
    for f in full_paths:
        reads_with_full_path = set()
        for read_id, path in tree.find_all(f):
            reads_with_full_path.add(read_id.replace("_reverse", ""))
        if len(reads_with_full_path) >= threshold:
            full_blocks[tuple(f)] = reads_with_full_path


def filter_blocks(full_blocks):
    filtered_blocks = {}
    for p in sorted(list(full_blocks.keys()), key=len, reverse=True):
        p_list = list(p)
        rv_p_list = list(reversed(p_list))
        if not any(
            is_sublist(list(f), p_list) or is_sublist(list(f), rv_p_list) for f in filtered_blocks
        ):
            filtered_blocks[p] = full_blocks[p]
    return filtered_blocks


def get_unique_anchor_suffixes(tree, anchor, node_mapping, anchor_suffixes):
    suffixes = get_suffixes_from_initial_tree(tree, anchor)
    for read_id in suffixes:
        path_list = suffixes[read_id]
        AMR_indices = [1 if n in node_mapping else 0 for n in path_list]
        if len(path_list) > 1:
            if AMR_indices[1] == 1:
                if anchor not in anchor_suffixes:
                    anchor_suffixes[anchor] = {}
                path_tuple = tuple(path_list)
                if path_tuple not in anchor_suffixes[anchor]:
                    anchor_suffixes[anchor][path_tuple] = set()
                anchor_suffixes[anchor][path_tuple].add(read_id.replace("_reverse", ""))


def filter_anchor_suffixes(anchor_suffixes, threshold):
    filtered_anchor_suffixes = {}
    for a in anchor_suffixes:
        for f in sorted(list(anchor_suffixes[a].keys()), key=len, reverse=True):
            if len(anchor_suffixes[a][f]) >= threshold:
                if a not in filtered_anchor_suffixes:
                    filtered_anchor_suffixes[a] = {}
                filtered_anchor_suffixes[a][f] = anchor_suffixes[a][f]
    return filtered_anchor_suffixes
