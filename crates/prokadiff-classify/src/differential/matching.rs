use std::cmp::Ordering;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct CoordinateCost {
    pub(super) max_delta: i64,
    pub(super) sum_delta: i64,
}

impl CoordinateCost {
    const ZERO: Self = Self {
        max_delta: 0,
        sum_delta: 0,
    };

    fn checked_add(self, other: Self) -> Option<Self> {
        Some(Self {
            max_delta: self.max_delta.checked_add(other.max_delta)?,
            sum_delta: self.sum_delta.checked_add(other.sum_delta)?,
        })
    }

    const fn reversed(self) -> Self {
        Self {
            max_delta: -self.max_delta,
            sum_delta: -self.sum_delta,
        }
    }
}

impl Ord for CoordinateCost {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.max_delta, self.sum_delta).cmp(&(other.max_delta, other.sum_delta))
    }
}

impl PartialOrd for CoordinateCost {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct TolerantCandidate {
    pub(super) edited_index: usize,
    pub(super) starter_index: usize,
    pub(super) cost: CoordinateCost,
}

#[derive(Clone, Copy, Debug)]
struct ResidualEdge {
    to: usize,
    reverse_index: usize,
    capacity: u8,
    cost: CoordinateCost,
}

fn add_edge(graph: &mut [Vec<ResidualEdge>], from: usize, to: usize, cost: CoordinateCost) {
    let forward_index = graph[from].len();
    let reverse_index = graph[to].len();
    graph[from].push(ResidualEdge {
        to,
        reverse_index,
        capacity: 1,
        cost,
    });
    graph[to].push(ResidualEdge {
        to: from,
        reverse_index: forward_index,
        capacity: 0,
        cost: cost.reversed(),
    });
}

pub(super) fn minimum_cost_maximum_matching(
    edited_count: usize,
    starter_count: usize,
    candidates: &[TolerantCandidate],
) -> Vec<(usize, usize)> {
    let source = 0;
    let edited_start = 1;
    let starter_start = edited_start + edited_count;
    let sink = starter_start + starter_count;
    let mut graph = vec![Vec::new(); sink + 1];

    for edited_index in 0..edited_count {
        add_edge(
            &mut graph,
            source,
            edited_start + edited_index,
            CoordinateCost::ZERO,
        );
    }
    for candidate in candidates {
        add_edge(
            &mut graph,
            edited_start + candidate.edited_index,
            starter_start + candidate.starter_index,
            candidate.cost,
        );
    }
    for starter_index in 0..starter_count {
        add_edge(
            &mut graph,
            starter_start + starter_index,
            sink,
            CoordinateCost::ZERO,
        );
    }

    loop {
        let mut distance = vec![None; graph.len()];
        let mut previous = vec![None; graph.len()];
        distance[source] = Some(CoordinateCost::ZERO);

        for _ in 1..graph.len() {
            let mut changed = false;
            for node in 0..graph.len() {
                let Some(prefix_cost) = distance[node] else {
                    continue;
                };
                for (edge_index, edge) in graph[node].iter().enumerate() {
                    if edge.capacity == 0 {
                        continue;
                    }
                    let Some(candidate_cost) = prefix_cost.checked_add(edge.cost) else {
                        continue;
                    };
                    if distance[edge.to].is_none_or(|current| candidate_cost < current) {
                        distance[edge.to] = Some(candidate_cost);
                        previous[edge.to] = Some((node, edge_index));
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
        }

        if distance[sink].is_none() {
            break;
        }

        let mut node = sink;
        while node != source {
            let Some((from, edge_index)) = previous[node] else {
                return collect_matches(&graph, edited_start, starter_start, starter_count);
            };
            let reverse_index = graph[from][edge_index].reverse_index;
            graph[from][edge_index].capacity = 0;
            graph[node][reverse_index].capacity = 1;
            node = from;
        }
    }

    collect_matches(&graph, edited_start, starter_start, starter_count)
}

fn collect_matches(
    graph: &[Vec<ResidualEdge>],
    edited_start: usize,
    starter_start: usize,
    starter_count: usize,
) -> Vec<(usize, usize)> {
    let starter_end = starter_start + starter_count;
    let mut matches = Vec::new();
    for (edited_index, edges) in graph[edited_start..starter_start].iter().enumerate() {
        for edge in edges {
            if (starter_start..starter_end).contains(&edge.to) && edge.capacity == 0 {
                matches.push((edited_index, edge.to - starter_start));
            }
        }
    }
    matches
}
