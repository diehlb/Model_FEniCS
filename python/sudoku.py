#!/opt/local/bin/python
#  Originally an assignment for an udacity class on software testing the following code generates
#  a possible solution to a sudoku puzzle if one exists.
#
#  Currently there is NO nice method of inputting the puzzle, as seen at the end in the testing portion 
#  the input is a list of 9 lists of 9 integers where 0 is a blank.
#  
#  To approach this problem a class Puzzle was created, a fuction 'check_sudoku' 
#(from previous assignment in the class) to check the validity of a current puzzle,
# and then a function 'solve_sudoku' that picks squares to guess possible solutions = brute force
#
import time

def check_sudoku(grid):
# grid is a list of len 9 with each element being a row (length 9) of the sudoku puzzle
# Three possible returns
#       -None   The input 'grid' is not of valid type, list of 9 rows each having 9 elements 0-9 
#       -False  The non-zero entries in the grid contain an obvious inconsistency
#       -True   The non-zero entries seem consistent
    try:
        if type(grid) != list or len(grid) != 9: return None
    except:
        return None

    for row in grid:
        try:
            if type(row) != list or len(row) != 9: return None
        except:
            return None
        for e in row:
            if type(e) != int or e < 0 or e > 9: return None        #Invalid cell
        #Catch duplicates in row
        chk = chkNonet(row)
        if not chk: return chk
    
    for i in range(9):
        #Catch duplicates in col
        col = [row[0] for row in grid]          
        chk = chkNonet(col)
        if not chk: return chk

        #Catch duplicates in nonet
        #   A nonet being one of the 3x3 subsets in the puzzle
        nonet = [x for row in grid[3*(i/3):3*(i/3)+3] for x in row[3*(i%3):3*(i%3)+3]]
        chk = chkNonet(nonet)
        if not chk: return chk
    
    return True

class Puzzle:
#  A class is created for a sudoku puzzle to track the possibilities for given cells as well as
#  current entries
#  The possible entries are 1-9 following game rules, with 0 indicating an empty cell
    def __init__(self, data):
        self.grid     = data
        self.avail_row = [ [1]*10 for i in range(9)]   # boolean possibilities per row
        self.avail_col = [ [1]*10 for i in range(9)]   #    "          "       per col
        self.avail_non = [ [1]*10 for i in range(9)]   #    "          "       per nonet
        self.nAvail   = [9]*27                        # rows, cols, then nonets

# On intitialization the given grid is examined to fill in the other data, mainly 'available' numbers
        for i in range(9):
            row = data[i]
            col = [R[i] for R in data]
            non = [x for R in data[3*(i/3):3*(i/3)+3] for x in R[3*(i%3):3*(i%3)+3]]
            for j in range(9):
                # for any number in the grid, it cannot be in any other spot in the row, col, or nonet
                self.avail_row[i][row[j]] = 0
                self.avail_col[i][col[j]] = 0
                self.avail_non[i][non[j]] = 0
            self.countAvail()       #This counts the possibilities to determine where guesses will occur

    def inc(self, ind):
    #Using brute force, the cell being guessed at is incremented if previous guess was proven inconsistent
    #This is more than just adding one, as it will check that availability in the row, col and nonet
        r = ind/9  
        c = ind%9
        n = (c/3) + 3*(r/3)
        old = self.grid[r][c]
        if old > 0:
            self.zero(ind)
        new = 10                #Initialize in way to catch null possibilities
        
        for i in range(old+1,10):      
            if self.avail_row[r][i] & self.avail_col[c][i] & self.avail_non[n][i]:
                new = i
                break
        if new > 9:
            return False

        self.grid[r][c]= new

        self.avail_row[r][new] = 0
        self.avail_col[c][new] = 0
        self.avail_non[n][new] = 0

        return True

    def zero(self, ind):
    #If the brute force method found no possible consistent entry, 
    #zero it back out and back up to last cell
        r = ind/9
        c = ind%9
        n = (c/3) + 3*(r/3)
        old = self.grid[r][c]
        self.grid[r][c]= 0
        self.avail_row[r][old] = 1
        self.avail_col[c][old] = 1
        self.avail_non[n][old] = 1

    def getZ(self,grp):
    #The number of possibilities is a list joining the rows, cols, and nonets
    #When determining best cell to guess next it is a simple min of that large list,
    #so this function backs out the index (0-80) of the first zero in that group
        if  grp < 9:
            Nonet = self.grid[grp]
            Ind_non = Nonet.index(min(Nonet))
            return Ind_non + 9*grp
        elif grp < 18:
            Nonet = [row[grp-9] for row in self.grid]
            Ind_non = Nonet.index(min(Nonet))
            return (grp - 9) + 9*Ind_non
        else:
            i = grp - 18
            Nonet = [x for row in self.grid[3*(i/3):3*(i/3)+3] for x in row[3*(i%3):3*(i%3)+3]]
            Ind_non = Nonet.index(min(Nonet))
            return Ind_non%3 + 3*(i%3) + 9*(Ind_non/3) + 27*(i/3)

    def countAvail(self):
    #Adjust the tally for possibilities of the rows, cols, and nonets
        for i in range(9):
            self.nAvail[i] = self.avail_row[i][1:].count(1)
            self.nAvail[i+9] = self.avail_col[i][1:].count(1)
            self.nAvail[i+18] = self.avail_non[i][1:].count(1)
        for j in range(27):
            #When the row,col,nonet is full it has zero available
            #But to prevent looking when using min, the value is bumped to 10
            # (above the 9 possible ints)
            if self.nAvail[j] == 0:
               self.nAvail[j] = 10

#####

def chkNonet(nonet):
#Given a list this creates a dictionary to check no nonzero value is repeated
    d = {}
    chk = True
    for val in nonet:
        if val != 0 and val in d:
            chk = False; return chk
        d[val] = 1
    return chk

def printGrid(grid):
#A simple print function so solutions are displayed in a readable way
    for row in grid:
        print row

def solve_sudoku(grid):
    if not check_sudoku(grid):  #Check sudoku grid is valid
        return None
    puz   = Puzzle(grid)        #Initialize 
    fills = []                  #Create fill list

    #nAvail = 10 means the group is filled, so this while checks that
    #there is at least one nonfilled group: else it's done
    while (min(puz.nAvail) < 10):
         grp = puz.nAvail.index(min(puz.nAvail))    #Find the most complete group
         cel = puz.getZ(grp)                        #Get the index of a cell in the group
         assert puz.grid[cel/9][cel%9] == 0         #  Check it is blank
         assert max(puz.nAvail) < 11                #  Check our possibilities count is sane
          
         fills.append(cel)
         valid = puz.inc(cel)                       #  Increment the cell

         while (not valid) and len(fills) > 1:      #  If cell can't be incremented and there is a prev
             cel = fills.pop()                      #guess, then back up 
             puz.zero(cel)
             valid = puz.inc(fills[-1])

         if len(fills) == 1 and (not valid):        #  If cell can't be incremented and there is NO prev
             return False                           #there is no solution

         puz.countAvail()                           # Update possibilities

    assert check_sudoku(puz.grid)                   # Be careful and check each guess is consistent

    return puz.grid 
       


##########################
# The following is a bit of code to check the solver and time some modifications
# First some test cases (taken from an udacity class assignment) followed by a function and script
evil = [[6,3,0,0,0,2,0,0,5 ],
        [0,0,2,9,0,0,0,1,0 ],
        [0,0,9,0,0,6,0,0,0 ],
        [0,0,5,8,0,0,1,0,0 ],
        [3,0,0,0,1,0,0,0,9 ],
        [0,0,6,0,0,7,5,0,0 ],
        [0,0,0,4,0,0,9,0,0 ],
        [0,8,0,0,0,5,2,0,0 ],
        [4,0,0,6,0,0,0,0,7 ]]

ill_formed = [[5,3,4,6,7,8,9,1,2],
              [6,7,2,1,9,5,3,4,8],
              [1,9,8,3,4,2,5,6,7],
              [8,5,9,7,6,1,4,2,3],
              [4,2,6,8,5,3,7,9],  # <---
              [7,1,3,9,2,4,8,5,6],
              [9,6,1,5,3,7,2,8,4],
              [2,8,7,4,1,9,6,3,5],
              [3,4,5,2,8,6,1,7,9]]

# check_sudoku should return True
Valid = [[5,3,4,6,7,8,9,1,2],
         [6,7,2,1,9,5,3,4,8],
         [1,9,8,3,4,2,5,6,7],
         [0,0,0,7,6,1,0,0,3],
         [0,0,0,8,5,3,0,0,1],
         [0,0,0,9,2,4,0,0,6],
         [9,6,1,5,3,7,2,8,4],
         [2,8,7,4,1,9,6,3,5],
         [3,4,5,2,8,6,1,7,9]]

# check_sudoku should return False
invalid = [[5,3,4,6,7,8,9,1,2],
           [6,7,2,1,9,5,3,4,8],
           [1,9,8,3,8,2,5,6,7],
           [8,5,9,7,6,1,4,2,3],
           [4,2,6,8,5,3,7,9,1],
           [7,1,3,9,2,4,8,5,6],
           [9,6,1,5,3,7,2,8,4],
           [2,8,7,4,1,9,6,3,5],
           [3,4,5,2,8,6,1,7,9]]

# check_sudoku should return True
easy = [[2,9,0,0,0,0,0,7,0],
        [3,0,6,0,0,8,4,0,0],
        [8,0,0,0,4,0,0,0,2],
        [0,2,0,0,3,1,0,0,7],
        [0,0,0,0,8,0,0,0,0],
        [1,0,0,9,5,0,0,6,0],
        [7,0,0,0,9,0,0,0,1],
        [0,0,1,2,0,0,3,0,6],
        [0,3,0,0,0,0,0,5,9]]

# check_sudoku should return True
hard = [[1,0,0,0,0,7,0,9,0],
        [0,3,0,0,2,0,0,0,8],
        [0,0,9,6,0,0,5,0,0],
        [0,0,5,3,0,0,9,0,0],
        [0,1,0,0,8,0,0,0,2],
        [6,0,0,0,0,4,0,0,0],
        [3,0,0,0,0,0,0,1,0],
        [0,4,0,0,0,0,0,0,7],
        [0,0,7,0,0,0,3,0,0]]



def test_check():
    t0 = time.time()
    print solve_sudoku(ill_formed) # --> None
    printGrid(solve_sudoku(Valid))     # --> True
    print solve_sudoku(invalid)    # --> False
    printGrid(solve_sudoku(hard))       # --> True
    print '-----------------------------'
    printGrid(solve_sudoku(easy))       # --> True
    t1 = time.time()
    T = (t1-t0) * 1000
    print 'Tests took %0.4f ms' % T

t0 = time.time()
res = solve_sudoku(hard)
if type(res) == list:
    printGrid(res)
else: print res
t1 = time.time()
T = (t1 - t0) * 1000
print 'Test took %0.4f ms' % T
