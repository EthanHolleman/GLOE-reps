# Update the workflow readme image
import os

if not os.path.isdir('resources'):
    os.mkdir('resources')

image_path = 'resources/current_workflow.png'

print('Generating DAG image')
status = os.system(
    f'snakemake -j 1 --dag | dot -Tpng > {image_path}'
)

# read README.md and replace current workflow image, if not there
# append it to the end
print('Updating README.md')
with open('README.md', 'r+') as readme:
    content = readme.read()
    readme.seek(0)
    if os.path.exists(image_path) and image_path not in content:
        readme.write(
            f'\n# Current Workflow\n![]({image_path})'
        )
        os.system(f'git add README.md {image_path}')
        os.system('git commit -m "update workflow image"')
        os.system('git push')
print('Done')


