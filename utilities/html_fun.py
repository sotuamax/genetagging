html_content = """
<!DOCTYPE html>
<html>
<head>
    <title>My Page</title>
</head>
<body>
    <h1>Hello from Python!</h1>
    <p>This is a paragraph.</p>
</body>
</html>
"""

# Write to a file
with open("output.html", "w") as file:
    file.write(html_content)

    