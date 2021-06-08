print('The top of the greeting_module has been read.')

def hello(name):
  greeting = "Hello {}!".format(name)
  return greeting
 
def ahoy(name):
  greeting = "Ahoy {}!".format(name)
  return greeting

x = 5

print('The bottom of the greeting_module has been read.')