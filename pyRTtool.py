import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import matplotlib as mpl

mpl.rc('font',size=17)


class RT:
    def __init__(self, prefix="model", directory="sims"):
        self.prefix = prefix
        self.directory = directory
        self.coolants = ['CO', 'C+', 'C', 'O']  # Default coolants
        self.additional_coolants = []
        self.all_coolants = self.coolants.copy()
        
        # Initialize transition map for each coolant
        self.transition_map = {
            'CO': OrderedDict([(f'CO{i+1}{i}', i) for i in range(10)]),
            'C+': OrderedDict([('C+10', 0)]),
            'C': OrderedDict([('CI10', 0), ('CI21', 1)]),
            'O': OrderedDict([('OI63', 0), ('OI146', 1)]),
            'HCO+': OrderedDict([(f'HCO+{i+1}{i}', i) for i in range(10)])
        }
        
        # Check for additional coolants by looking for files in the directory
        self._detect_additional_coolants()
        
    def _detect_additional_coolants(self):
        """Detect additional coolants by scanning the directory for files"""
        if not os.path.exists(self.directory):
            return
            
        files = os.listdir(self.directory)
        for f in files:
            if f.startswith(f"Tr_{self.prefix}.") and f.endswith(".dat"):
                coolant = f.split('.')[1]
                if coolant not in self.coolants and coolant not in self.additional_coolants:
                    self.additional_coolants.append(coolant)
                    # Initialize transition map for new coolant (up to 10 transitions)
                    self.transition_map[coolant] = OrderedDict([(f'{coolant}{i+1}{i}', i) for i in range(10)])
        
        self.all_coolants.extend(self.additional_coolants)
    
    def _read_file(self, filename):
        """Read all lines of a file"""
        try:
            with open(os.path.join(self.directory, filename), 'r') as f:
                return [line.strip() for line in f.readlines() if line.strip()]
        except FileNotFoundError:
            return None
    
    def show(self, coolant):
        """Display the last row of data for a specific coolant"""
        if coolant not in self.all_coolants:
            print(f"Coolant '{coolant}' not found. Available coolants: {', '.join(self.all_coolants)}")
            return
            
        # Read Tr file
        tr_file = f"Tr_{self.prefix}.{coolant}.dat"
        tr_lines = self._read_file(tr_file)
        
        # Read tau file
        tau_file = f"tau_{self.prefix}.{coolant}.dat"
        tau_lines = self._read_file(tau_file)
        
        if tr_lines is None or tau_lines is None:
            print(f"Could not find data files for coolant '{coolant}'")
            return
            
        # Get the last line of each file
        tr_line = tr_lines[-1]
        tau_line = tau_lines[-1]
        
        # Parse the data
        tr_parts = tr_line.split()
        tau_parts = tau_line.split()
        
        # Get the number of transitions
        num_transitions = len(self.transition_map.get(coolant, {}))
        
        # Extract the values
        Ntot = tr_parts[0]
        NH2 = tr_parts[1]
        Tgas = tr_parts[2]
        Ncoolant = tr_parts[3]
        Tr_values = tr_parts[4:4+num_transitions]
        tau_values = tau_parts[:num_transitions]
        
        # Print the results
        print(f"N({coolant})= {Ncoolant} cm-2")
        print(" --------------------")
        print("      Line            tau         Tr (K)")
        
        # Get transition labels in correct order (1-0, 2-1, ..., 10-9)
        transitions = list(self.transition_map.get(coolant, {}).keys())
        
        for i, trans in enumerate(transitions):
            if i < len(Tr_values) and i < len(tau_values):
                # Format the transition name based on coolant type
                if coolant == 'CO':
                    # Format as "CO (J-J')" where J is upper level
                    upper_level = int(trans[2:-1])
                    formatted_trans = f"CO ({upper_level}-{upper_level-1})"
                elif coolant == 'C+':
                    # Only one transition: [CII] 158um
                    formatted_trans = "[CII] 158um"
                elif coolant == 'C':
                    # Format as [CI] (J-J')
                    if trans == 'CI10':
                        formatted_trans = "[CI] (1-0)"
                    else:
                        formatted_trans = "[CI] (2-1)"
                elif coolant == 'O':
                    # Format as [OI] wavelength
                    if trans == 'OI63':
                        formatted_trans = "[OI] 63um"
                    else:
                        formatted_trans = "[OI] 146um"
                elif coolant == 'HCO+':
                    # Format as "HCO+ (J-J')" where J is upper level
                    # Extract the numbers after HCO+ and before the last character
                    numbers = trans[4:-1]  # Gets the part between 'HCO+' and the last character
                    if numbers.isdigit():
                        upper_level = int(numbers)
                        formatted_trans = f"HCO+ ({upper_level}-{upper_level-1})"
                    else:
                        formatted_trans = trans  # Fallback to original if parsing fails
                else:
                    # For additional coolants, use the raw transition name
                    formatted_trans = trans
                
                print(f"{formatted_trans:>15} :  {tau_values[i]:>12}  {Tr_values[i]:>12}")


    def plot(self, x_var, y_var, scale='loglog'):
        """
        Plot column density vs Tr or tau.
        
        Parameters:
        x_var: str - Column density variable ('Ntot', 'NH2', or 'Ncoolant' like 'NCO', 'NC+', etc.)
        y_var: str - Transition variable ('Tr(CO10)', 'tau(HCO+10)', 'Tr(C+)', etc.)
        scale: str - Scale type ('loglog', 'semilogx', 'semilogy', 'linear')
        """
        # Parse y_var to get coolant and transition
        if y_var.startswith('Tr('):
            y_type = 'Tr'
            trans_str = y_var[3:-1]
        elif y_var.startswith('tau('):
            y_type = 'tau'
            trans_str = y_var[4:-1]
        else:
            raise ValueError("y_var must start with 'Tr(' or 'tau('")
        
        # Special case for C+ which only has one transition
        if trans_str == 'C+':
            coolant = 'C+'
            transition = 0
        else:
            # Parse the transition string to get coolant and transition
            coolant = None
            transition = None
            
            # Check for known coolant patterns
            for c in self.transition_map:
                if trans_str.startswith(c):
                    coolant = c
                    transition = self.transition_map[c].get(trans_str)
                    break
        
        if coolant is None or transition is None:
            raise ValueError(f"Could not parse transition '{trans_str}'. For C+, use 'Tr(C+)' or 'tau(C+)'")
        
        # Read the data files
        tr_file = f"Tr_{self.prefix}.{coolant}.dat"
        tr_lines = self._read_file(tr_file)
        
        if y_type == 'tau':
            tau_file = f"tau_{self.prefix}.{coolant}.dat"
            tau_lines = self._read_file(tau_file)
            if tau_lines is None:
                raise FileNotFoundError(f"Could not find tau file for coolant '{coolant}'")
        
        if tr_lines is None:
            raise FileNotFoundError(f"Could not find Tr file for coolant '{coolant}'")
        
        # Prepare data for plotting
        x_data = []
        y_data = []
        
        for line in tr_lines:
            parts = line.split()
            # Get x value
            if x_var == 'Ntot':
                x_val = float(parts[0])
            elif x_var == 'NH2':
                x_val = float(parts[1])
            elif x_var.startswith('N') and x_var[1:] in self.all_coolants:
                if parts[3] == 'NaN':
                    continue  # Skip NaN values
                x_val = float(parts[3])
            else:
                raise ValueError(f"Unrecognized x_var '{x_var}'. Use 'Ntot', 'NH2', or 'Ncoolant' like 'NCO', 'NC+', etc.")
            
            # Get y value
            if y_type == 'Tr':
                y_val = float(parts[4 + transition])
            else:
                # For tau, we need to find the corresponding line in the tau file
                idx = tr_lines.index(line)
                if idx < len(tau_lines):
                    tau_parts = tau_lines[idx].split()
                    if transition < len(tau_parts):
                        y_val = float(tau_parts[transition])
                    else:
                        continue  # Skip if transition not available
                else:
                    continue  # Skip if no corresponding tau line
            
            x_data.append(x_val)
            y_data.append(y_val)
        
        # Create the plot
        plt.figure(figsize=(8, 6))
        
        if scale == 'loglog':
            plt.loglog(x_data, y_data, 'o-')
        elif scale == 'semilogx':
            plt.semilogx(x_data, y_data, 'o-')
        elif scale == 'semilogy':
            plt.semilogy(x_data, y_data, 'o-')
        elif scale == 'linear':
            plt.plot(x_data, y_data, 'o-')
        else:
            raise ValueError(f"Unrecognized scale '{scale}'. Use 'loglog', 'semilogx', 'semilogy', or 'linear'")
        
        # Set labels
        x_label_map = {
            'Ntot': 'Total column density (cm-2)',
            'NH2': 'H2 column density (cm-2)',
        }
        
        if x_var.startswith('N') and x_var[1:] in self.all_coolants:
            x_label = f'{x_var[1:]} column density (cm-2)'
        else:
            x_label = x_label_map.get(x_var, x_var)
        
        y_label = f"{y_var.split('(')[0]}({trans_str})"
        
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(f"{x_label} vs {y_label}")
        plt.grid(True, which="both", ls="-")
        plt.tight_layout()
        plt.show()
    
    def list_coolants(self):
        """List all available coolants"""
        print("Available coolants:")
        for coolant in self.all_coolants:
            print(f" - {coolant}")
    
    def list_transitions(self, coolant):
        """List available transitions for a coolant"""
        if coolant not in self.all_coolants:
            print(f"Coolant '{coolant}' not found. Available coolants: {', '.join(self.all_coolants)}")
            return
        
        transitions = self.transition_map.get(coolant, {})
        if not transitions:
            print(f"No transition information available for '{coolant}'")
            return
        
        print(f"Available transitions for {coolant}:")
        for trans in transitions.keys():
            print(f" - {trans}")

    def sled(self, coolant, norm=None, scale='log'):
        """
        Plot Spectral Line Energy Distribution (SLED) for a coolant.

        Parameters:
        -----------
        coolant : str
            Coolant name (e.g., 'CO', 'C+', 'HCO+')
        norm : str or None, optional
            Transition to normalize to (e.g., 'CO10', 'HCO+21'). If None, no normalization.
        scale : str, optional
            Y-axis scale ('linear' or 'log'). Default is 'log'.
        """
        if coolant not in self.all_coolants:
            print(f"Coolant '{coolant}' not found. Available coolants: {', '.join(self.all_coolants)}")
            return

        # Read Tr file
        tr_file = f"Tr_{self.prefix}.{coolant}.dat"
        tr_lines = self._read_file(tr_file)

        if tr_lines is None:
            print(f"Could not find Tr file for coolant '{coolant}'")
            return

        # Get the last row data
        last_line = tr_lines[-1].split()
        Ncoolant = last_line[3]
        Tr_values = [float(x) for x in last_line[4:4+len(self.transition_map.get(coolant, {}))]]

        # Prepare x-axis labels and positions
        x_labels = []
        x_positions = []

        # Special cases for C+ and O
        if coolant == 'C+':
            x_labels = ['158μm']
            x_positions = [1]
        elif coolant == 'O':
            x_labels = ['63μm', '146μm']
            x_positions = [1, 2]
        elif coolant == 'C':
            x_labels = ['[CI](1-0)', '[CI](2-1)']
            x_positions = [1, 2]
        else:
            # For CO, HCO+, and other coolants
            transitions = sorted(self.transition_map.get(coolant, {}).items(), key=lambda x: x[1])
            x_positions = [i+1 for i in range(len(transitions))]
            if coolant == 'CO':
                x_labels = [f'CO({i+1}-{i})' for i in range(len(transitions))]
            else:
                x_labels = [f'{coolant}({i+1}-{i})' for i in range(len(transitions))]

        # Apply normalization if requested
        if norm is not None:
            # Find the index of the transition to normalize to
            norm_index = None
            for i, (trans, idx) in enumerate(self.transition_map.get(coolant, {}).items()):
                if trans == norm:
                    norm_index = idx
                    break

            if norm_index is not None and norm_index < len(Tr_values):
                norm_value = Tr_values[norm_index]
                if norm_value != 0:
                    Tr_values = [val/norm_value for val in Tr_values]
                    ylabel = f'Tr/Tr({norm})'
                else:
                    print(f"Cannot normalize to {norm} (zero value)")
            else:
                print(f"Normalization transition '{norm}' not found for {coolant}")
        else:
            ylabel = 'Tr (K)'

        # Create the plot
        plt.figure(figsize=(10, 6))

        # Plot the points
        plt.plot(x_positions, Tr_values, 'o-', markersize=8, linewidth=2)

        # Set scales
        if scale == 'log':
            plt.yscale('log')
        elif scale != 'linear':
            print(f"Warning: Unknown scale '{scale}'. Using 'linear'")

        # Customize the plot
        plt.xticks(x_positions, x_labels, rotation=45, ha='right')
        plt.xlabel('Transition')
        plt.ylabel(ylabel)

        title = f"{coolant} Spectral Line Energy Distribution"
        if norm is not None:
            title += f" (normalized to {norm})"
        plt.title(title)

        plt.grid(True, which='both', linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.show()

    def makeRT(self):
        """
        Run 'make clean; make' in the RT-tool/ directory.
        """
        original_dir = os.getcwd()
        os.chdir("RT-tool")
        os.system("make clean; make")
        os.chdir(original_dir)
    
    def runRT(self):
        """
        Execute ./RTtool in the current directory.
        """
        os.system('./RTtool')
    
    def params(self, directory, prefix, vturb=1.0):
        """
        Create a paramsRT.dat file with the specified parameters.
        
        Parameters:
        directory (str): Directory path (mandatory)
        prefix (str): File prefix (mandatory)
        vturb (float): Turbulent velocity (default=1.0)
        """
        try:
            with open("paramsRT.dat", "w") as f:
                f.write(f"{directory}\n")
                f.write(f"{prefix}\n")
                f.write(f"{vturb}\n")
            
            print("Successfully created paramsRT.dat")
            return True
            
        except Exception as e:
            print(f"Error creating paramsRT.dat: {str(e)}")
            return False
