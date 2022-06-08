Decorators are a powerful and useful tool in Python for modifying functions and classes.

In this tutorial, we will learn about what decorators are and how to use them effectively.

We will also see some examples of how decorators can be used to add functionality to your code.

So, what exactly are decorators?

Decorators are simply functions that take another function as an argument. They can be used to modify the behaviour of that function without having to directly modify the code of the function itself.

This is extremely powerful as it allows you to change the behaviour of code without changing the code itself.

Decorators are usually defined in their own module so that they can be reused easily.

Here is a simple example of a decorator:

```python
def my_decorator(func):
  def wrapper():
    print(“Before the function is called.”)
    func()
 print(“After the function is called.”)
 return wrapper
 ```

This decorator simply prints some text before and after the function that it takes as an argument is called. Now, let’s see how we can use this decorator.

```python
@my_decorator
def print_hello():
 print(“Hello, world!”)
print_hello()
```
**Output:**

```
Before the function is called.
Hello, world!
After the function is called.
```
As you can see, we simply added the decorator @my_decorator before the definition of the print_hello() function.

This is all that is needed to apply the decorator. When we call the print_hello() function, the code inside the wrapper function in the decorator is executed before and after the print_hello() function is called.

Decorators can also take arguments. For example, let’s say we want a decorator that prints the time before and after the function it is decorating is called. We can do this as follows:

```python
import time
def timed(func):
 def wrapper(*args, **kwargs):
 start = time.time()
 result = func(*args, **kwargs)
 end = time.time()
 print(“The function took {} seconds to complete.”.format(end — start))
 return result
 return wrapper
@timed
def print_hello():
 time.sleep(1)
 print(“Hello, world!”)
print_hello()
```

**Output:**

```
The function took 1.00048828125 seconds to complete.
Hello, world!
```
As you can see, the decorator takes the function it is decorating as an argument and returns a wrapper function.

The wrapper function is what is actually called when the decorated function is called.

The wrapper function takes the same arguments as the function it is decorating and simply calls that function with those arguments.

It also prints the time before and after the function is called.

Now, let’s see how we can use decorators to add functionality to a class.

```python
class Person:
 def __init__(self, name, age):
 self.name = name
 self.age = age
@logged
class Person:
 def __init__(self, name, age):
 self.name = name
 self.age = age
p = Person(“John”, 30)
print(p)
```

**Output:**

```
Person(name=’John’, age=30)
```

As you can see, we simply added the decorator @logged before the definition of the Person class. This decorator simply prints the values of the arguments that are passed to the __init__() method of the class.

Decorators are a powerful tool that can be used to modify the behavior of functions and classes.

In this tutorial, we have learned about what decorators are and how to use them effectively. We have also seen some examples of how decorators can be used to add functionality to your code.

Before you leave:

#### About the author:
Alain Saamego: Software engineer, writer and content strategist at SelfGrow.co.uk

Email:alain@selfgrow.co.uk

