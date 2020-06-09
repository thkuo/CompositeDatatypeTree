
def outer(func):
    def multi(*args):
        return('{}.{}'.format(args[0], args[0]))
    def dec(*args):
        return('{}-WHEEE-{}'.format(args[0], args[0]))
    def opt(*args):
        if func(*args) == 'a':
#            print('This is type a')
            return(multi(*args))
        elif func(*args) == 'b':
#            print('This is type b')
            return(dec(*args))
        else:
            return('what...')
    return(opt)


@outer
def determine_func(x):
    d= ['a', 'b', 'c']
    return(d[x%len(d)])

print(determine_func(0))
print(determine_func(1))
print(determine_func(2))
