# Sudoku Generator
## Usage
```pip install numpy matplotlyb```  
```python sudoku.py```  
then Sudoku pazzle and an solution will be generated as *pazzle.png*, *answer.png*
## Logic
Define a function which samples multiple solutions from a given hints. (hints are the pre-locked numbers) 
Using those solutions, make a new hint. 
Starting from zero hints, execute that function multiple times until generating only one solution. 
Enjoy the pazzle 🥸
