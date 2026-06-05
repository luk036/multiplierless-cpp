import sys

with open(sys.argv[1], "r") as f:
    lines = f.readlines()

fixed = []
in_port_list = False

for line in lines:
    s = line.strip()
    if "module" in s and "(" in s:
        in_port_list = True
        fixed.append(line)
        continue
    if in_port_list and ");" in s:
        in_port_list = False
        if fixed and fixed[-1].rstrip().endswith(","):
            fixed[-1] = fixed[-1].rstrip().rstrip(",") + "\n"
        fixed.append(line)
        continue
    if in_port_list and s.startswith("output") and not s.endswith(","):
        fixed.append(line.rstrip() + ",\n")
    else:
        fixed.append(line)

with open(sys.argv[1], "w") as f:
    f.writelines(fixed)

print(f"Fixed {sys.argv[1]}")
