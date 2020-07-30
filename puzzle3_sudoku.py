from sympy.utilities.iterables import multiset_permutations
import numpy as np

class GiantSudoku:
    def __init__(self):
        self.matrix = np.zeros((9,9))
        self.x = [None]*9
        self.y = [None]*9
        self.dx = np.array([[1,2,125,1,29,1,1,8,1],
                            [22,1,1,1,79,1,1,2,1]])
        self.dy = np.array([[4,1,2,1,1,11,127,1,4],
                            [257,1,17,1,1,5,877,1,1]])

        self.perms = self.getPermutations()
        self.initSearchSpaces()
        self.reduce()

    def getPermutations(self):
        permGen = multiset_permutations(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]))
        perms = np.array([perm for perm in permGen])
        pot10 = np.power(10, np.arange(8,-1,-1))
        nums = np.sum(pot10 * perms, axis=-1)
        return nums

    def initSearchSpaces(self):
        for i in range(9):
            D1, D2 = self.dx[:,i]
            self.x[i] = self.getDivisiblePermutations(D1, D2)
            D1, D2 = self.dy[:, i]
            self.y[i] = self.getDivisiblePermutations(D1, D2)

    def getDivisiblePermutations(self, D1, D2):
        perms = np.copy(self.perms)
        perms = perms[perms % D1 == 0]
        restrictedPerms = []
        if D2 == 0: return np.copy(perms)
        for num in perms:
            if self.invertInt(num) % D2 == 0:
                restrictedPerms.append(num)
        return np.array(restrictedPerms)

    def invertInt(self, num):
        return int(str(num)[::-1])

    def reduceMagicSquare(self):
        # center of magic sqr is 5
        self.matrix[4,4] = 5
        for row in range(3,6):
            maskX = []
            for i in range(len(self.x[row])):
                sum = self.intAt(self.x[row][i], 3)
                sum+= self.intAt(self.x[row][i], 4)
                sum+= self.intAt(self.x[row][i], 5)
                maskX.append(sum==15)
            self.x[row] = self.x[row][np.array(maskX)]

        for col in range(3, 6):
            maskY = []
            for i in range(len(self.y[col])):
                sum = self.intAt(self.y[col][i], 3)
                sum += self.intAt(self.y[col][i], 4)
                sum += self.intAt(self.y[col][i], 5)
                maskY.append(sum == 15)
            self.y[col] = self.y[col][np.array(maskY)]

    def intAt(self, num, pos):
        return num//10**(8-pos) % 10

    def reduceColRow(self):
        for i in range(9):
            self.x[i] = self.reduceRow(i)
            self.y[i] = self.reduceCol(i)

    def reduceCol(self, col):
        new = []
        for entry in self.y[col]:
            valid = True
            cause = None
            for pos in range(9):
                d = self.intAt(entry, pos)
                dmap = self.matrix[pos, col]
                if dmap != 0 and dmap != d:
                    valid = False
                if d in self.matrix[pos, col + 1:] or d in self.matrix[pos, :col]:
                    valid = False
                boxContainsEntry = self.isAtBox(pos, col, d)
                if boxContainsEntry:
                    valid = False
            if valid:
                new.append(entry)
        return new

    def reduceRow(self, row):
        new = []
        for entry in self.x[row]:
            valid = True
            cause = None
            for pos in range(9):
                d = self.intAt(entry, pos)
                dmap = self.matrix[row, pos]
                if dmap != 0 and dmap != d:
                    valid = False
                if d in self.matrix[row + 1:, pos] or d in self.matrix[:row, pos]:
                    valid = False
                boxContainsEntry = self.isAtBox(row, pos, d)
                if boxContainsEntry:
                    valid = False
            if valid:
                new.append(entry)
        return new

    def isAtBox(self, row, col, entry):
        infX = 3 * (row // 3)
        infY = 3 * (col // 3)
        relXPos = row % 3
        relYPos = col % 3

        box = np.copy(self.matrix[infX:infX + 3, infY:infY + 3])
        box[relXPos, relYPos] = 0
        if entry in box:
            return True
        return False

    def reduceUnmatch(self, row, col):
        xval = []
        for entry in self.x[row]:
            xval.append(self.intAt(entry, col))
        xval = set(xval)

        yval = []
        for entry in self.y[col]:
            yval.append(self.intAt(entry, row))
        yval = set(yval)

        xincomp = xval - yval
        yincomp = yval - xval

        newx = []
        for entry in self.x[row]:
            d = self.intAt(entry, col)
            if d not in xincomp:
                newx.append(entry)

        newy = []
        for entry in self.y[col]:
            d = self.intAt(entry, row)
            if d not in yincomp:
                newy.append(entry)

        return newx, newy

    def reduceUnmatches(self):
        for row in range(9):
            for col in range(9):
                if len(self.x[row]) > 0 and len(self.y[col]) > 0:
                    self.x[row], self.y[col] = self.reduceUnmatch(row, col)


    def markInvariants(self):
        for col in range(9):
            for row in range(9):
                if self.isInvariant(self.y[col], row):
                    if len(self.y[col]) > 0:
                        if self.matrix[row, col] == 0:
                            inv = self.intAt(self.y[col][0], row)
                            self.matrix[row, col] = inv
                    else:
                        self.matrix[row, col] = self.intAt(self.y[col][0],row)
                if self.isInvariant(self.x[row], col):
                    if len(self.x[row]) > 1:
                        if self.matrix[row, col] == 0:
                            inv = self.intAt(self.x[row][0], col)
                            self.matrix[row, col] = inv
                    else:
                        self.matrix[row, col] = self.intAt(self.x[row][0],col)

    def isInvariant(self, group, i):
        freq = []
        for entry in group:
            freq.append(self.intAt(entry, i))
        freq = set(freq)
        return (len(freq) == 1)

    def countMarked(self):
        return np.count_nonzero(self.matrix)

    def getSpaceSize(self):
        sz = 0
        for i in range(9):
            sz += len(self.x[i])
            sz += len(self.y[i])
        return sz

    def reduce(self):
        sz = self.getSpaceSize()
        mk = self.countMarked()
        self.reduceMagicSquare()
        while(mk<81):
            lastSz = sz
            lastMk = mk
            self.reduceUnmatches()
            self.reduceColRow()
            self.markInvariants()
            sz = self.getSpaceSize()
            mk = self.countMarked()
            print(f"marked {mk-lastMk} new cells. Total = {mk}/81")
            print(f"reduced space by {lastSz -sz}. From {lastSz} to {sz} with. min=18")
            if lastSz == sz:
                break


if __name__ == "__main__":
    sdk = Sudoku()
    print(sdk.matrix)