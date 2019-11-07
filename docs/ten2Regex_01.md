#### Credit: this tutorial is extracted from [Regex tutorial — A quick cheatsheet by examples](https://medium.com/factory-mind/regex-tutorial-a-simple-cheatsheet-by-examples-649dc1c3f285)

Regular expressions (regex or regexp) are extremely useful in **extracting information from any text** by searching for one or more matches of a specific search pattern (i.e. a specific sequence of ASCII or unicode characters).

Fields of application range from validation to parsing/replacing strings, passing through translating data to other formats and web scraping.

One of the most interesting features is that once you’ve learned the syntax, you can actually use this tool in (almost) all programming languages ​​(JavaScript, Java, VB, C #, C / C++, Python, Perl, Ruby, Delphi, R, Tcl, and many others) with the slightest distinctions about the support of the most advanced features and syntax versions supported by the engines).

Let’s start by looking at some examples and explanations.


## Basic topics


### Anchors — ^ and $


<table>
  <tr>
   <td><strong><code>^The</code></strong>
   </td>
   <td>matches any string that <strong>starts with</strong> <strong>The </strong>-><strong><a href="https://regex101.com/r/cO8lqs/2"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><strong><code>end$</code></strong>
   </td>
   <td>matches <strong>a</strong> string that <strong>ends with</strong> <strong>end</strong>
   </td>
  </tr>
  <tr>
   <td><strong><code>^The end$</code></strong>
   </td>
   <td><strong>exact string match</strong> (starts and <strong>ends</strong> with <strong>The end</strong>)
   </td>
  </tr>
  <tr>
   <td><strong><code>roar</code></strong>
   </td>
   <td>matches any string that <strong>has the text roar in it</strong>
   </td>
  </tr>
</table>



### Quantifiers — * + ? and {}


<table>
  <tr>
   <td><code>abc<strong>*</strong></code>
   </td>
   <td>matches a string that has <strong>ab followed by zero or more c </strong>-><strong><a href="https://regex101.com/r/cO8lqs/1"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><code>abc<strong>+</strong></code>
   </td>
   <td>matches a string that has <strong>ab followed by one or more c</strong>
   </td>
  </tr>
  <tr>
   <td><code>abc<strong>?</strong></code>
   </td>
   <td>matches a string that has <strong>ab followed by zero or one c</strong>
   </td>
  </tr>
  <tr>
   <td><code>abc<strong>{2}</strong></code>
   </td>
   <td>matches a string that has <strong>ab followed by 2 c</strong>
   </td>
  </tr>
  <tr>
   <td><code>abc<strong>{2,}</strong></code>
   </td>
   <td>matches a string that has <strong>ab followed by 2 or more c</strong>
   </td>
  </tr>
  <tr>
   <td><code>abc<strong>{2,5}</strong></code>
   </td>
   <td>matches a string that has <strong>ab followed by 2 up to 5 c</strong>
   </td>
  </tr>
  <tr>
   <td><code>a<strong>(bc)*</strong></code>
   </td>
   <td>matches a string that has <strong>a followed by zero or more copies of the sequence bc</strong>
   </td>
  </tr>
  <tr>
   <td><code>a<strong>(bc){2,5}</strong></code>
   </td>
   <td>matches a string that has<strong> a followed by 2 up to 5 copies of the sequence bc</strong>
   </td>
  </tr>
</table>



### OR operator — | or []


<table>
  <tr>
   <td><strong><code>a(b|c)</code></strong>
   </td>
   <td>matches a string that has <strong>a followed by b or c -><a href="https://regex101.com/r/cO8lqs/3"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><strong><code>a[bc]</code></strong>
   </td>
   <td>same as previous
   </td>
  </tr>
</table>



### ` `Character classes — \d \w \s and .


<table>
  <tr>
   <td><strong><code>\d</code></strong>
   </td>
   <td>matches a <strong>single character</strong> that is a <strong>digit </strong>-><strong><a href="https://regex101.com/r/cO8lqs/4"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><strong><code>\w</code></strong>
   </td>
   <td>matches a <strong>word character</strong> (alphanumeric character plus underscore) -><a href="https://regex101.com/r/cO8lqs/4"> Try it!</a>
   </td>
  </tr>
  <tr>
   <td><strong><code>\s</code></strong>
   </td>
   <td>matches a <strong>whitespace character</strong> (includes tabs and line breaks)
   </td>
  </tr>
  <tr>
   <td><strong><code>.</code></strong>
   </td>
   <td>matches <strong>any character </strong>-><strong><a href="https://regex101.com/r/cO8lqs/5"> Try it!</a></strong>
   </td>
  </tr>
</table>


Use the `.` operator carefully since often class or negated character class (which we’ll cover next) are faster and more precise.

`\d`, `\w` and `\s` also present their negations with `\D`, `\W` and `\S` respectively.

For example, `\D` will perform the inverse match with respect to that obtained with `\d`.


<table>
  <tr>
   <td><strong><code>\D</code></strong>
   </td>
   <td>matches a <strong>single non-digit character </strong>-><strong><a href="https://regex101.com/r/cO8lqs/6"> Try it!</a></strong>
   </td>
  </tr>
</table>


In order to be taken literally, you must escape the characters `^.[$()|*+?{\`with a backslash `\` as they have special meaning.


<table>
  <tr>
   <td><strong><code>\$\d</code></strong>
   </td>
   <td>matches a string that has a <strong>$ before one digit </strong>-><strong><a href="https://regex101.com/r/cO8lqs/9"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><strong><code>\$\d</code></strong>
   </td>
   <td>Notice that you can match also <strong>non-printable characters</strong> like tabs <code>\t</code>, new-lines <code>\n</code>, carriage returns <code>\r</code>.
   </td>
  </tr>
</table>



### Flags

We are learning how to construct a regex but forgetting a fundamental concept: **flags**.

A regex usually comes within this form **<code>/abc/</code></strong>, where the search pattern is delimited by two slash characters <code>/</code>. At the end we can specify a flag with these values (we can also combine them each other):



*   
**g **(global) does not return after the first match, restarting the subsequent searches from the end of the previous match


*   
**m** (multi-line) when enabled `^` and `$` will match the start and end of a line, instead of the whole string


*   
**i** (insensitive) makes the whole expression case-insensitive (for instance **<code>/aBc/i</code></strong> would match <strong><code>AbC</code></strong>)


---



## Intermediate topics


### Grouping and capturing — ()


<table>
  <tr>
   <td><code>a<strong>(</strong>bc<strong>)</strong></code>
   </td>
   <td>parentheses create a <strong>capturing group with value bc </strong>-><strong><a href="https://regex101.com/r/cO8lqs/11"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><code>a<strong>(?:</strong>bc<strong>)</strong>*</code>
   </td>
   <td>using <strong>?:</strong> we <strong>disable the capturing group </strong>-><strong><a href="https://regex101.com/r/cO8lqs/12"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><code>a<strong>(?<foo></strong>bc<strong>)</strong></code>
   </td>
   <td>using <strong>?<foo></strong> we put a name to the group -><a href="https://regex101.com/r/cO8lqs/17"> Try it!</a>
   </td>
  </tr>
</table>


This operator is very useful when we need to extract information from strings or data using your preferred programming language. Any multiple occurrences captured by several groups will be exposed in the form of a classical array: we will access their values specifying using an index on the result of the match.

If we choose to put a name to the groups (using <code>(<strong>?<foo></strong>...)</code>) we will be able to retrieve the group values using the match result like a dictionary where the keys will be the name of each group.


### Bracket expressions — []


<table>
  <tr>
   <td><strong><code>[abc]</code></strong>
   </td>
   <td>matches a string that has <strong>either an a or a b or a c </strong>-><strong> </strong>is the <strong>same as a|b|c </strong>-><strong><a href="https://regex101.com/r/cO8lqs/7"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><strong><code>[a-c]</code></strong>
   </td>
   <td>same as previous
   </td>
  </tr>
  <tr>
   <td><strong><code>[a-fA-F0-9]</code></strong>
   </td>
   <td>a string that represents <strong>a single hexadecimal digit, case insensitively </strong>-><strong><a href="https://regex101.com/r/cO8lqs/22"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><strong><code>[0-9]%</code></strong>
   </td>
   <td>a string that has a character <strong>from 0 to 9 before a % sign</strong>
   </td>
  </tr>
  <tr>
   <td><strong><code>[^a-zA-Z]</code></strong>
   </td>
   <td>a string that has <strong>not a letter from a to z or from A to Z. </strong>In this case the <strong>^</strong> is used as <strong>negation of the expression </strong>-><strong><a href="https://regex101.com/r/cO8lqs/10"> Try it!</a></strong>
   </td>
  </tr>
</table>



```


```


Remember that inside bracket expressions all special characters (including the backslash `\`) lose their special powers: thus we will not apply the “escape rule”.


### Greedy and Lazy match

The quantifiers ( `* + {}`) are greedy operators, so they expand the match as far as they can through the provided text.

For example, `<.+>` matches `<div>simple div</div>` in <code>This is a <strong><div> simple div</div></strong> test</code>. In order to catch only the <code>div</code> tag we can use a <code>?</code> to make it lazy:


<table>
  <tr>
   <td><strong><code><.+?></code></strong>
   </td>
   <td>matches <strong>any character one or more</strong> times included <strong>inside <</strong> and <strong>></strong>, <strong>expanding as needed </strong>-><strong><a href="https://regex101.com/r/cO8lqs/24"> Try it!</a></strong>
   </td>
  </tr>
</table>


Notice that a better solution should avoid the usage of `.` in favor of a more strict regex:


<table>
  <tr>
   <td><strong><code><[^<>]+></code></strong>
   </td>
   <td>matches <strong>any character except < or > one or more </strong>times included <strong>inside <</strong> and <strong>> </strong>-><strong><a href="https://regex101.com/r/cO8lqs/23"> Try it!</a></strong>
   </td>
  </tr>
</table>




---



## Advanced topics


### Boundaries — \b and \B


<table>
  <tr>
   <td><strong><code>\babc\b</code></strong>
   </td>
   <td>performs a <strong>"whole words only" search </strong>-><strong><a href="https://regex101.com/r/cO8lqs/25"> Try it!</a></strong>
   </td>
  </tr>
</table>


`\b` represents an **anchor like caret** (it is similar to `$` and `^`) matching positions where **one side is a word** **character **(like `\w`) and the **other side is not a word** **character **(for instance it may be the beginning of the string or a space character).

It comes with its **negation**, `\B`. This matches all positions where `\b` doesn’t match and could be if we want to find a search pattern fully surrounded by word characters.


<table>
  <tr>
   <td><strong><code>\Babc\B </code></strong>
   </td>
   <td>matches only if the pattern is <strong>fully surrounded by word</strong> characters -><a href="https://regex101.com/r/cO8lqs/26"> Try it!</a>
   </td>
  </tr>
</table>



```


```



### Back-references — \1


<table>
  <tr>
   <td><strong><code>([abc])\1</code></strong>
   </td>
   <td>using <strong>\1</strong> it matches <strong>the same</strong> text <strong>that was matched by the first capturing group </strong>-><strong><a href="https://regex101.com/r/cO8lqs/14"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><code>([abc])<strong>(</strong>[de]<strong>)\2</strong>\1</code>
   </td>
   <td>we can use <strong>\2</strong> (\3, \4, etc.) to identify <strong>the same</strong> text that <strong>was matched by the second </strong>(third, fourth, etc.) <strong>capturing group </strong>-><strong><a href="https://regex101.com/r/cO8lqs/15"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><strong><code>(?<foo>[abc])\k<foo></code></strong>
   </td>
   <td>we put the name <strong>foo t</strong>o the group and we reference it later (<strong>\k<foo></strong>). The result is the same of the first regex -><a href="https://regex101.com/r/cO8lqs/16"> Try it!</a>
   </td>
  </tr>
</table>



### Look-ahead and Look-behind — (?=) and (?<=)


<table>
  <tr>
   <td><code>d<strong>(?=</strong>r<strong>)</strong></code>
   </td>
   <td>matches a <strong>d </strong>only if is<strong> followed by r</strong>, <strong>but r will not be</strong> part of the overall regex <strong>match </strong>-><strong><a href="https://regex101.com/r/cO8lqs/18"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><strong><code>(?<=r)d</code></strong>
   </td>
   <td>matches a <strong>d </strong>only if is<strong> preceded by an r</strong>, <strong>but r will not be</strong> part of the overall regex <strong>match </strong>-><strong><a href="https://regex101.com/r/cO8lqs/19"> Try it!</a></strong>
   </td>
  </tr>
</table>


You can use also the negation operator!


<table>
  <tr>
   <td><code>d<strong>(?!</strong>r<strong>)</strong></code>
   </td>
   <td>matches a <strong>d </strong>only if is<strong> not followed by r</strong>, <strong>but r will not be</strong> part of the overall regex <strong>match </strong>-><strong><a href="https://regex101.com/r/cO8lqs/20"> Try it!</a></strong>
   </td>
  </tr>
  <tr>
   <td><strong><code>(?<!r)d</code></strong>
   </td>
   <td>matches a <strong>d </strong>only if is<strong> not preceded by an r</strong>, <strong>but r will not be</strong> part of the overall regex <strong>match </strong>-><strong><a href="https://regex101.com/r/cO8lqs/21"> Try it!</a></strong>
   </td>
  </tr>
</table>



```


```



## Summary

As you’ve seen, the application fields of regex can be multiple and I’m sure that you’ve recognized at least one of these tasks among those seen in your developer career, here a quick list:



*   
data validation (for example check if a time string i well-formed)


*   
data scraping (especially web scraping, find all pages that contain a certain set of words eventually in a specific order)


*   
data wrangling (transform data from “raw” to another format)


*   
string parsing (for example catch all URL GET parameters, capture text inside a set of parenthesis)


*   
string replacement (for example, even during a code session using a common IDE to translate a Java or C# class in the respective JSON object — replace “;” with “,” make it lowercase, avoid type declaration, etc.)


*   
syntax highlightning, file renaming, packet sniffing and many other applications involving strings (where data need not be textual)
Have fun and do not forget to recommend the article if you liked it 💚


<!-- Docs to Markdown version 1.0β17 -->