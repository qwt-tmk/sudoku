import copy
import numpy as np
import matplotlib.pyplot as plt


class Sudoku:
    def __init__(self, H: dict):
        self.rng = np.random.default_rng()
        self.fullset = {i for i in range(1, 10)}
        self.H = H  # { (i, j): number } <--- hint
        self.field = self.make_initial_field()
        self.solutions = []

    def make_initial_field(self):
        ar = np.zeros((9, 9), dtype=int)
        for p, v in self.H.items():
            ar[p] = v
        for i in range(9):
            s = set(ar[i, :]) - {0}
            comp = np.array(list(self.fullset - s))
            self.rng.shuffle(comp)
            k = 0
            for j in range(9):
                if ar[i, j] == 0:
                    ar[i, j] = comp[k]
                    k += 1
        return ar

    def get_block(self, field, k):
        r = 3 * (k // 3)
        c = 3 * (k % 3)
        return field[r : r + 3, c : c + 3].flatten().tolist()

    def get_column(self, field, j):
        return field[:, j].tolist()

    def calc_E_seg(self, seg):
        seg = set(seg)
        return len(self.fullset - seg)

    def calc_E(self, field):
        return sum(
            [
                self.calc_E_seg(self.get_column(field, j))
                + self.calc_E_seg(self.get_block(field, j))
                for j in range(9)
            ]
        )

    def calc_f(self, field, beta):
        return np.exp(-beta * self.calc_E(field))

    def get_two_points(self):
        i, j, j_ = self.rng.choice(range(9), 3)
        if (i, j) == (i, j_) or (i, j) in self.H or (i, j_) in self.H:
            return self.get_two_points()
        else:
            return (i, j), (i, j_)

    def swap(self, beta):
        p, p_ = self.get_two_points()
        canditate = self.field.copy()
        canditate[p], canditate[p_] = canditate[p_], canditate[p]
        diff = self.calc_E(canditate) - self.calc_E(self.field)
        if diff < 0:
            self.field = canditate.copy()
        else:
            if self.rng.random() < np.exp(-beta * diff):
                self.field = canditate.copy()
            else:
                pass

    def check(self):
        if self.calc_E(self.field) < 1:
            self.solutions.append(self.field.copy())

    def sample_sa(self, num_reads, MCS):
        for read in range(num_reads):
            self.field = self.make_initial_field()
            initial = 0.1
            end = 3.0
            beta_schedule_root = [
                (end - initial) * np.sqrt(m / MCS) + initial for m in range(MCS)
            ]
            for beta in beta_schedule_root:
                print(f"\r reads: {read+1}, beta: {beta}", end="")
                for t in range(1000):
                    self.swap(beta)
                    self.check()


def reduce_solutions(solutions: list) -> list:
    sol = copy.deepcopy(solutions)
    sol = np.array([arr.flatten() for arr in sol])
    _, indices = np.unique(sol, axis=0, return_index=True)
    return [solutions[i] for i in indices]


def make_hint(solutions: list[np.ndarray]) -> tuple:
    if len(solutions) < 2:
        raise ValueError(
            f"No needs to add a hint. number of solutions is {len(solutions)}"
        )
    accum = {((i, j), n): 0 for i in range(9) for j in range(9) for n in range(1, 10)}
    for sol in solutions:
        for idx, num in np.ndenumerate(sol):
            accum[(idx, int(num))] += 1

    mini = min(set(accum.values()) - {0})
    canditates = [idx for idx, v in accum.items() if v == mini]
    i = np.random.default_rng().choice(range(len(canditates)))
    return canditates[i]


def generate_pazzle(num_reads=100, MCS=1000):
    H = {}
    sol = np.zeros((9, 9)) - 1
    while True:
        print(f"num hints : {len(H)}")
        sudoku = Sudoku(H)
        print("sampling starts")
        sudoku.sample_sa(num_reads=num_reads, MCS=MCS)
        print("sampling ends")
        solutions = reduce_solutions(sudoku.solutions)
        num_sol = len(solutions)
        print(f"num solutions : {num_sol}")
        if num_sol == 1:
            sol = sudoku.solutions[0]
            break
        elif num_sol == 0:
            print("no solution found, pop one hint")
            _, _ = H.popitem()
        new_hint = make_hint(solutions)
        print(f"new hint : {new_hint}")
        H.update({new_hint[0]: new_hint[1]})
    return H, sol


def save_grid_fig(H: dict, file_path):
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot()
    for idx, num in H.items():
        ax.text(idx[0] + 0.5, idx[1] + 0.5, str(num), va="center", ha="center")
    # ax.set_xlim(range(10))
    # ax.set_ylim(range(10))
    ax.set_xticks(range(10))
    ax.set_yticks(range(10))
    ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
    ax.grid()
    fig.savefig(file_path)


if __name__ == "__main__":
    H, sol = generate_pazzle(num_reads=10, MCS=1000)
    dsol = {}
    for idx, num in np.ndenumerate(sol):
        dsol[idx] = num
    save_grid_fig(H, "pazzle.png")
    save_grid_fig(dsol, "answer.png")
