import cgen as c

func = c.FunctionBody(
    c.FunctionDeclaration(c.Const(c.Pointer(c.Value("char", "greet"))), []),
    c.Block([
        c.Statement('return "hello world"')
        ])
    )
code = c.Module([])

code.append(c.Value('int', 'cont'))
code.append(c.Assign('cont', '0'))
code.append(c.Increment('cont', '5'))


print(code)
