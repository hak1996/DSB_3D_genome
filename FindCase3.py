
def merge_point_pairs(pairs):
    merged_pairs = []
    while pairs:
        current_pair = pairs.pop(0)
        merged = False
        for pair in list(merged_pairs):
            if current_pair[0] in pair or current_pair[1] in pair:
                new_pair = tuple(sorted(set(current_pair + pair)))
                merged_pairs.remove(pair)
                merged_pairs.append(new_pair)
                merged = True
                break
        if not merged:
            merged_pairs.append(current_pair)
    return merged_pairs


class AdjacencyDatabase:
    def __init__(self):
        # 使用字典来存储邻接关系，键是点的id，值是与该点邻接的点的集合
        self.adjacency_dict = {}

    def add(self, point_id, adjacent_point_id):
        # 添加邻接关系
        if point_id not in self.adjacency_dict:
            self.adjacency_dict[point_id] = set()
        self.adjacency_dict[point_id].add(adjacent_point_id)

        # 如果是无向图，也需要添加反向的邻接关系
        if adjacent_point_id not in self.adjacency_dict:
            self.adjacency_dict[adjacent_point_id] = set()
        self.adjacency_dict[adjacent_point_id].add(point_id)

    def find_adjacent_pairs_in_points(self, points_to_check):
        # 在给定的待查点中找出所有满足邻接关系的点对
        adjacent_pairs = []
        for point_id in points_to_check:
            if point_id in self.adjacency_dict:
                for adjacent_point_id in self.adjacency_dict[point_id]:
                    if adjacent_point_id in points_to_check:
                        # 避免重复添加点对，只添加一次
                        if (adjacent_point_id, point_id) not in adjacent_pairs:
                            adjacent_pairs.append((point_id, adjacent_point_id))
        return adjacent_pairs

    def find_case3_cluster(self, points_to_check):
        adjacent_pairs = self.find_adjacent_pairs_in_points(points_to_check)
        merged_pairs = merge_point_pairs(adjacent_pairs)
        return merged_pairs





if __name__ == "__main__":
    # 使用示例
    db = AdjacencyDatabase()
    db.add(0, 1)
    db.add(1, 0)
    db.add(1, 2)  # 点1有多个邻接点
    db.add(2, 1)
    db.add(2,3)
    db.add(3,4)
    db.add(1,5)
    db.add(5,6)
    db.add(6,7)

    # 找出所有满足邻接关系的点对
    points = [0,2,3,4,1,6,7]
    adjacent_pairs = db.find_adjacent_pairs_in_points(points)
    print(adjacent_pairs)
    cluster = db.find_case3_cluster(points)
    print(cluster)
