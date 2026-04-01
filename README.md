# Minimum Cost Flow - Orlin (1993)

Υλοποιηση βασισμενη στο paper:
"A Faster Strongly Polynomial Minimum Cost Flow Algorithm"
του James B. Orlin, Operations Research, Vol. 41, No. 2 (1993)

Περιεχει δυο εκδοσεις:
- **Section 3**: The Edmonds-Karp Scaling Technique
- **Section 4**: The Strongly Polynomial Time Algorithm

## Compilation

```
make
```

## Εκτελεση

```
./mcf_section3 < input.txt
./mcf_section4 < input.txt
```

## Μορφη αρχειου εισοδου

```
n m
b(0) b(1) ... b(n-1)
from0 to0 cost0
from1 to1 cost1
...
```

- **Γραμμη 1**: `n` = αριθμος κομβων, `m` = αριθμος ακμων
- **Γραμμη 2**: Οι τιμες supply/demand καθε κομβου χωρισμενες με κενο (θετικο = supply, αρνητικο = demand, 0 = transit)
- **Γραμμες 3 εως m+2**: Καθε ακμη ως `from to cost`

### Παραδειγμα

```
4 5
10 0 0 -10
0 1 2
0 2 5
1 2 1
1 3 8
2 3 3
```

4 κομβοι, 5 ακμες. Ο κομβος 0 εχει supply 10, ο κομβος 3 εχει demand 10.

## Generator

Δημιουργια τυχαιων δικτυων δοκιμων:

```
python generator.py <nodes> <sources> <sinks> <supply> <output> [--seed N]
```

Παραδειγμα:
```
python generator.py 10 2 2 100 network.txt
```
