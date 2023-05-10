# Optimization Methods for Cybersecurity
# Campopiano Daniele - 174624

import cplex
from cplex.exceptions import CplexError
import os


# POPOLAMENTO DEL PROBLEMA PER COLONNA
def populate_by_column(my_prob, my_A, my_Aeq, my_c, my_lb, my_ub, my_ctype, my_colnames, my_b, my_beq, my_sense,
                       my_rownames, my_num_cols):
    for j in range(my_num_cols):
        variables = []
        coefficients = []
        lower_bounds = []
        upper_bounds = []
        types = []
        for i in range(len(my_A)):
            variables.append(i)
            coefficients.append(my_A[i][j])
            lower_bounds.append(0.0)
            upper_bounds.append(1.0)
            types.append(my_ctype[j])
        for i in range(len(my_Aeq)):
            variables.append(i + len(my_A))
            coefficients.append(my_Aeq[i][j])
            lower_bounds.append(0.0)
            upper_bounds.append(1.0)
            types.append(my_ctype[j])
        colname = "x_" + str(j)
        try:
            my_prob.variables.add(obj=my_c[j], lb=my_lb[j], ub=my_ub[j], types=my_ctype[j], names=[colname])
        except TypeError:
            print(f"Error: my_c[{j}]={my_c[j]} is not a valid input for CPLEX variables.add() function.")
            break
        my_prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=variables, val=coefficients)], senses=my_sense,
                                       rhs=[my_b], names=[colname])
        my_prob.quadratic_constraints.add(lin_expr=[cplex.SparsePair(ind=variables, val=coefficients)],
                                          quad_expr=[[], [], []], sense='E', rhs=[0.0], names=[colname])


# POPOLAMENTO DEL PROBLEMA PER RIGA
def populate_by_row(my_prob, my_A, my_Aeq, my_c, my_lb, my_ub, my_ctype, my_colnames, my_b, my_beq, my_sense,
                    my_rownames):
    if len(my_A) != len(my_b) or len(my_A[0]) != len(my_c):
        print("Error: dimension mismatch between A, b and c")
        return
    i = 0
    rows = []
    cols = []
    vals = []

    numcols = len(my_c)
    numrows = len(my_b)

    for i in range(numrows):
        row = []
        row_eq = []
        for j in range(numcols):
            val = my_A[i][j]
            if val != 0:
                rows.append(i)
                cols.append(j)
                vals.append(val)
                row.append(j)
        if my_sense[i] == 'E':
            my_prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=row, val=vals)],
                                           senses=[my_sense[i]], rhs=[my_beq[i]], names=[my_rownames[i]])
            print("1")
        else:
            my_prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=row, val=vals)],
                                           senses=[my_sense[i]], rhs=[my_b[i]], names=[my_rownames[i]])
            print("2")
        vals.clear()
        row.clear()

    for i in range(len(my_Aeq)):
        row_eq = []
        for j in range(numcols):
            val = my_Aeq[i][j]
            if val != 0:
                row_eq.append(j)
        my_prob.quadratic_constraints.add(lin_expr=[cplex.SparsePair(ind=row_eq, val=my_Aeq[i])],
                                          rhs=[my_beq[i]], sense='E', name='qc' + str(i + 1))
        row_eq.clear()

    for j in range(numcols):
        if my_ctype[j] == 'C':
            my_prob.variables.add(obj=[my_c[j]], lb=[my_lb[j]], ub=[my_ub[j]], types=[my_ctype[j]],
                                  names=[my_colnames[j]])
        else:
            my_prob.variables.add(obj=[my_c[j]], lb=[my_lb[j]], ub=[my_ub[j]], names=[my_colnames[j]])


# POPOLAMENTO DEL PROBLEMA PER COEFFIIENTI DIVERSI DA ZERO
def populate_by_nonzero(my_prob, my_A, my_Aeq, my_c, my_lb, my_ub, my_ctype, my_colnames, my_b, my_beq, my_sense, my_rownames):
    # Verifica che my_beq non sia vuota
    if not my_beq:
        return

    # Verifica che la lunghezza di my_beq sia maggiore o uguale a numrows
    numrows = len(my_sense)
    if len(my_beq) < numrows:
        return

    # Popola il problema con coefficienti diversi da zero
    for i in range(numrows):
        if my_sense[i] == 'E':
            my_prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=my_Aeq[i], val=[1.0]*len(my_Aeq[i]))],
                                            rhs=my_beq[i], names=[my_rownames[i]])
        else:
            my_prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=my_A[i], val=[1.0]*len(my_A[i]))],
                                            rhs=my_b[i], names=[my_rownames[i]])
    for j in range(len(my_c)):
        if my_c[j] != 0.0:
            my_prob.objective.set_linear([(my_colnames[j], my_c[j])])

#
def resolution_mps(pop_method, my_prob):
    try:
        # Estrae rhs, il vettore beq e b
        beq = my_prob.quadratic_constraints.get_rhs()  # Estrae i termini quadrati
        b = my_prob.linear_constraints.get_rhs()  # Estrae i termini lineari
        c = my_prob.objective.get_linear()  # Estrae la funzione obiettivo

        # Estrae i limiti inferiori e superiori delle variabili
        lb = my_prob.variables.get_lower_bounds()  # Estrae i limiti inferiori
        ub = my_prob.variables.get_upper_bounds()  # Estrae i limiti superiori

        # Estrae le informazioni di colonne e righe
        numcols = my_prob.variables.get_num()  # Estrae il numero di colonne
        colnames = my_prob.variables.get_names()  # Estrae i nomi delle colonne
        numrows = my_prob.linear_constraints.get_num()  # Estrae il numero di righe
        rownames = my_prob.linear_constraints.get_names()  # Estrae i nomi delle righe

        # Estrae ctype e sense
        ctype = my_prob.variables.get_types()  # Estrae il tipo di variabile
        sense = my_prob.linear_constraints.get_senses()  # Estrae il tipo di vincolo

        rows = []
        cols = []
        vals = []
        A = []
        Aeq = []

        # Estrae la matrice del vincolo in forma di triple (riga, colonna, valore)
        for i in range(numrows):
            row = my_prob.linear_constraints.get_rows(i)
            for j in range(len(row.ind)):
                rows.append(i)
                cols.append(row.ind[j])
                vals.append(row.val[j])
            if sense[i] == 'E':
                Aeq.append(vals.copy())
            else:
                A.append(vals.copy())
            vals.clear()

        # Popola il problema di ottimizzazione a seconda del metodo scelto
        if pop_method == "r":
            populate_by_row(my_prob, A, Aeq, c, lb, ub, ctype, colnames, b, beq, sense, rownames)
        elif pop_method == "c":
            populate_by_column(my_prob, A, Aeq, c, lb, ub, ctype, colnames, b, beq, sense, rownames, numcols)
        elif pop_method == "n":
            populate_by_nonzero(my_prob, A, Aeq, c, lb, ub, ctype, colnames, b, beq, sense, rownames)
        else:
            raise ValueError('pop_method must be one of "r", "c" or "n"')

        my_prob.solve()  # Risolve il problema di ottimizzazione

    except CplexError as exc:
        print(exc)
        return

    print("Solution status = ", my_prob.solution.get_status(), ":", end=' ')
    print(my_prob.solution.status[my_prob.solution.get_status()])
    print("Solution value  = ", my_prob.solution.get_objective_value())

    slack = my_prob.solution.get_linear_slacks()
    x = my_prob.solution.get_values()

    for j in range(numrows):
        print("Row %d:  Slack = %10f" % (j, slack[j]))
    for j in range(numcols):
        print("Column %d:  Value = %10f" % (j, x[j]))


if __name__ == "__main__":

    myProb = cplex.Cplex()

    mps_file = "Mps Files/blend2.mps"

    if os.path.exists(mps_file):
        print("\nAnalizzo il file nella path " + mps_file)
    else:
        print("\nFile non trovato nella path " + mps_file)

    # COLONNA
    myProb.read(mps_file)
    print("\n      c          generate problem by column\n")
    resolution_mps("c", myProb)

    # RIGA
    myProb.read(mps_file)
    print("\n      r          generate problem by row\n")
    resolution_mps("r", myProb)

    # COEFFICIENTI DIVERSI DA ZERO
    myProb.read(mps_file)
    print("\n      n          generate problem by nonzero\n")
    resolution_mps("n", myProb)
