while True:
    
    expression = input('In: ')
    
    if expression == '0':
        break
    elif expression == '':
        continue

    if(expression[-1] == ';'):
        expression = expression[:-1]

    if 'for ' in expression:
        for i in range(len(expression)):
            if expression[i] == '(':
                expression = expression[i+1:-1]
                expression = expression.split('; ')
                print('c.For(\'%s\', \'%s\', \'%s\', \nc.Block([\n\n])\n)' %(expression[0], expression[1], expression[2]))
                break
    elif expression[:3] == 'int' or expression[:4] == 'void' or expression[:6] == 'double':
        
        type_ = expression.split(' ')[0]
        expression = expression.strip(type_+' ')
        
        if(' = ' in expression):
            expression_split = expression.split(' = ')
            print('c.Initializer(c.Value(\'%s\', \'%s\'), \'%s\')' %(type_, expression_split[0], expression_split[1]))
        else:
            l_params = []
            for i in range(len(expression)):
                if expression[i] == '(':
                    name = expression[:i]
                    params = expression[i+1:-1]
                    params = params.split(', ')
                    for j in params:
                        for k in j.split(' '):
                            l_params.append(k)
                    break
            params_gen = ''
            for el in range(0, len(l_params), 2):
                params_gen += 'c.Value(\'%s\', \'%s\'), ' %(l_params[el], l_params[el+1])
            params_gen = params_gen[:-2]
            # print('type_ ', type_)
            # print('name ', name)
            # print('l_params ', l_params)
            # break
            print('code.append(c.FunctionBody(\nc.FunctionDeclaration(c.Value(%s, %s), [%s]), \nc.Block([\n\n])\n))' %(type_, name, params_gen))
    elif(' = ' in expression):
        expression_split = expression.split(' = ')
        print('c.Assign(\'%s\', \'%s\')' %(expression_split[0], expression_split[1]))